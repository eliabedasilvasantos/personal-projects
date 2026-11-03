"""Cold barotropic equations of state. Public dimensional values are SI.

rho0 is baryonic rest-mass density; epsilon includes rest-mass energy.
No built-in parameter set is claimed to fit nuclear matter or observations.
Tables are interpolated as log epsilon(log P) with monotone PCHIP and an
explicit, user-selected low-density extension. Density jumps are rejected.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from hashlib import sha256
import json
from math import asinh, exp, expm1, isfinite, log, log1p, pi, sqrt
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from .constants import C, HBAR, M_NEUTRON, MEV_FM3_TO_PA
from .exceptions import EOSDomainError


@dataclass(frozen=True)
class EOSState:
    pressure_pa: float
    energy_density_j_m3: float
    rest_density_kg_m3: float
    sound_speed_squared_over_c2: float


def _nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not isfinite(value) or value < 0:
        raise EOSDomainError(f"{name} must be nonnegative and finite")
    return value


class EquationOfState(ABC):
    """Protocol used by the integrator; enthalpy is dimensionless.

    state_from_h(0) returns the interior zero-pressure limit. Implementations
    must not silently extrapolate a finite table's high-density boundary.
    """

    name: str
    max_enthalpy: float = float("inf")
    validation_only: bool = False

    @abstractmethod
    def state_from_h(self, h: float) -> EOSState:
        pass

    @abstractmethod
    def h_from_rest_density(self, rho0: float) -> float:
        pass

    @abstractmethod
    def metadata(self) -> dict[str, Any]:
        pass

    def h_from_pressure(self, pressure_pa: float) -> float:
        p = _nonnegative(pressure_pa, "pressure")
        if p == 0:
            return 0.0
        hi = min(0.125, self.max_enthalpy)
        while self.state_from_h(hi).pressure_pa < p:
            if hi >= self.max_enthalpy or hi >= 30:
                raise EOSDomainError("Pressure exceeds the supported EOS domain")
            hi = min(2 * hi, self.max_enthalpy)
        # Relative comparison avoids an absolute root tolerance in Pa.
        return float(brentq(lambda h: self.state_from_h(h).pressure_pa / p - 1,
                            0, hi, xtol=5e-324, rtol=1e-13))

    def state_from_rest_density(self, rho0: float) -> EOSState:
        return self.state_from_h(self.h_from_rest_density(rho0))


@dataclass(frozen=True)
class RelativisticPolytrope(EquationOfState):
    """P=K rho0^Gamma; epsilon=rho0 c²+P/(Gamma-1).

    K has SI units Pa/(kg m^-3)^Gamma, which depend on Gamma. Prefer the
    from_reference constructor when comparing exponents at fixed density.
    """

    K: float
    gamma: float
    name: str = "Cold relativistic polytrope"

    def __post_init__(self) -> None:
        if not isfinite(self.K) or self.K <= 0:
            raise ValueError("K must be positive and finite")
        if not isfinite(self.gamma) or self.gamma <= 1:
            raise ValueError("gamma must be >1 and finite")

    @classmethod
    def from_reference(cls, rho_ref: float, pressure_ref: float, gamma: float,
                       name: str = "Reference-normalized polytrope") -> "RelativisticPolytrope":
        if min(rho_ref, pressure_ref) <= 0:
            raise ValueError("Reference density and pressure must be positive")
        return cls(pressure_ref / rho_ref**gamma, gamma, name)

    def h_from_rest_density(self, rho0: float) -> float:
        rho = _nonnegative(rho0, "rest density")
        return log1p(self.gamma * self.K * rho**(self.gamma-1) /
                     ((self.gamma-1)*C**2))

    def state_from_h(self, h: float) -> EOSState:
        h = _nonnegative(h, "enthalpy")
        if h == 0:
            return EOSState(0., 0., 0., 0.)
        rho = (expm1(h)*(self.gamma-1)*C**2/(self.gamma*self.K))**(1/(self.gamma-1))
        p = self.K * rho**self.gamma
        eps = rho*C*C + p/(self.gamma-1)
        return EOSState(p, eps, rho, self.gamma*p/(eps+p))

    def h_from_pressure(self, pressure_pa: float) -> float:
        p = _nonnegative(pressure_pa, "pressure")
        return self.h_from_rest_density((p/self.K)**(1/self.gamma))

    def metadata(self) -> dict:
        return {"type":"polytrope", "name":self.name, "K_SI":self.K,
                "gamma":self.gamma, "status":"uncalibrated phenomenological model"}


@dataclass(frozen=True)
class DegenerateNeutronGas(EquationOfState):
    """Exact zero-temperature free spin-1/2 neutron gas, degeneracy g=2.

    Nuclear interactions, leptons, beta equilibrium and a crust are absent.
    Low-x series prevent catastrophic cancellation in analytic Fermi integrals.
    """

    neutron_mass_kg: float = M_NEUTRON
    name: str = "Free degenerate neutron gas"

    def __post_init__(self) -> None:
        if not isfinite(self.neutron_mass_kg) or self.neutron_mass_kg <= 0:
            raise ValueError("Particle mass must be positive")

    @property
    def rho_scale(self) -> float:
        return self.neutron_mass_kg**4 * C**3 / (3*pi*pi*HBAR**3)

    def h_from_rest_density(self, rho0: float) -> float:
        x = (_nonnegative(rho0, "rest density")/self.rho_scale)**(1/3)
        return 0.5*log1p(x*x)

    def state_from_h(self, h: float) -> EOSState:
        h = _nonnegative(h, "enthalpy")
        x = sqrt(expm1(2*h))
        if x == 0:
            return EOSState(0., 0., 0., 0.)
        B = self.neutron_mass_kg**4 * C**5 / (pi*pi*HBAR**3)
        if x < 0.08:
            t = x*x
            # Integrals of x² sqrt(1+x²) and x⁴/sqrt(1+x²), through x^15.
            eps = B*x**3*(1/3+t*(1/10+t*(-1/56+t*(1/144+t*(-5/1408+
                       t*(7/3328-t*21/15360))))))
            p = B/3*x**5*(1/5+t*(-1/14+t*(1/24+t*(-5/176+
                     t*(35/1664-t*63/3840)))))
        else:
            root = sqrt(1+x*x)
            eps = B/8*(x*root*(1+2*x*x)-asinh(x))
            p = B/24*(x*root*(2*x*x-3)+3*asinh(x))
        return EOSState(p, eps, self.rho_scale*x**3, x*x/(3*(1+x*x)))

    def metadata(self) -> dict:
        return {"type":"fermi", "name":self.name, "neutron_mass_kg":self.neutron_mass_kg,
                "degeneracy":2, "temperature_K":0,
                "status":"noninteracting reference, not realistic neutron-star matter"}


class PiecewisePolytrope(EquationOfState):
    """Continuous cold piecewise polytrope, including energy offsets.

    Given K0, exponents and rest-density breaks, K_i ensure pressure
    continuity. epsilon=(1+a_i)rho0 c²+P/(Gamma_i-1) and a_i ensure energy
    and chemical-potential continuity. No density-jump phase transition.
    This implements the form, NOT any named Read et al. parameter fit.
    """

    def __init__(self, K0: float, gammas: list[float] | tuple[float, ...],
                 breaks_kg_m3: list[float] | tuple[float, ...],
                 name: str = "Continuous piecewise polytrope") -> None:
        self.name = name
        self.gammas = np.asarray(gammas, dtype=float)
        self.breaks = np.asarray(breaks_kg_m3, dtype=float)
        if len(self.gammas) != len(self.breaks)+1 or not len(self.gammas):
            raise ValueError("Need one more exponent than density breaks")
        if (not np.all(np.isfinite(self.gammas)) or np.any(self.gammas <= 1) or
            not np.all(np.isfinite(self.breaks)) or np.any(self.breaks <= 0) or
            np.any(np.diff(self.breaks) <= 0) or not isfinite(K0) or K0 <= 0):
            raise ValueError("Invalid piecewise-polytrope parameters")
        self.Ks = [float(K0)]
        self.offsets = [0.0]
        for i, rho in enumerate(self.breaks):
            old_g, new_g = self.gammas[i:i+2]
            new_k = self.Ks[-1]*rho**(old_g-new_g)
            p = self.Ks[-1]*rho**old_g
            a = self.offsets[-1]+p/(rho*C*C)*(1/(old_g-1)-1/(new_g-1))
            self.Ks.append(float(new_k))
            self.offsets.append(float(a))
        if any(1+a <= 0 for a in self.offsets):
            raise ValueError("This parametrization needs 1+a_i>0 in every segment")
        self.h_breaks = np.array([self._h_at_rho(rho, i) for i, rho in enumerate(self.breaks)])

    def _h_at_rho(self, rho: float, i: int) -> float:
        g, k, a = self.gammas[i], self.Ks[i], self.offsets[i]
        return log1p(a+g/(g-1)*k*rho**(g-1)/C**2)

    def h_from_rest_density(self, rho0: float) -> float:
        rho = _nonnegative(rho0, "rest density")
        return self._h_at_rho(rho, int(np.searchsorted(self.breaks, rho, side="right")))

    def state_from_h(self, h: float) -> EOSState:
        h = _nonnegative(h, "enthalpy")
        if h == 0:
            return EOSState(0., 0., 0., 0.)
        i = int(np.searchsorted(self.h_breaks, h, side="right"))
        g, k, a = self.gammas[i], self.Ks[i], self.offsets[i]
        rho = ((expm1(h)-a)*(g-1)*C*C/(g*k))**(1/(g-1))
        p = k*rho**g
        eps = (1+a)*rho*C*C+p/(g-1)
        return EOSState(float(p), float(eps), float(rho), float(g*p/(eps+p)))

    def metadata(self) -> dict:
        return {"type":"piecewise", "name":self.name, "K0_SI":self.Ks[0],
                "gammas":self.gammas.tolist(), "breaks_kg_m3":self.breaks.tolist(),
                "energy_offsets":self.offsets, "status":"uncalibrated, continuous, no density jumps"}


@dataclass(frozen=True)
class ConstantSoundSpeed(EquationOfState):
    """Self-bound linear EOS P=s(epsilon-epsilon_surface), 0<s<=1.

    rho0_surface is an explicit baryon-mass normalization. It does not follow
    from the linear P(epsilon) relation. No hadron/quark transition is implied.
    """

    s: float
    epsilon_surface_j_m3: float
    rho0_surface_kg_m3: float
    name: str = "Self-bound constant-sound-speed parametrization"

    def __post_init__(self) -> None:
        if not 0 < self.s <= 1 or not isfinite(self.s):
            raise ValueError("s = c_s²/c² must be in (0,1]")
        if (min(self.epsilon_surface_j_m3, self.rho0_surface_kg_m3) <= 0 or
            not isfinite(self.epsilon_surface_j_m3+self.rho0_surface_kg_m3)):
            raise ValueError("Surface energy and rest densities must be positive and finite")

    def h_from_rest_density(self, rho0: float) -> float:
        rho = _nonnegative(rho0, "rest density")
        if rho < self.rho0_surface_kg_m3:
            raise EOSDomainError("Self-bound material has no interior states below its surface density")
        return self.s*log(rho/self.rho0_surface_kg_m3)

    def h_from_pressure(self, pressure_pa: float) -> float:
        p = _nonnegative(pressure_pa, "pressure")
        return self.s/(1+self.s)*log1p((1+self.s)*p/(self.s*self.epsilon_surface_j_m3))

    def state_from_h(self, h: float) -> EOSState:
        h = _nonnegative(h, "enthalpy")
        p = self.epsilon_surface_j_m3*self.s/(1+self.s)*expm1((1+self.s)*h/self.s)
        eps = self.epsilon_surface_j_m3+p/self.s
        rho = self.rho0_surface_kg_m3*exp(h/self.s)
        return EOSState(p, eps, rho, self.s)

    def metadata(self) -> dict:
        return {"type":"css", "name":self.name, "s":self.s,
                "epsilon_surface_j_m3":self.epsilon_surface_j_m3,
                "rho0_surface_kg_m3":self.rho0_surface_kg_m3,
                "status":"uncalibrated self-bound surrogate, baryon normalization explicit"}


@dataclass(frozen=True)
class IncompressibleEOS(EquationOfState):
    """Acausal constant-energy-density exact-solution test ONLY.

    No physical baryon density is supplied: rho0=NaN disables Mb. Must call
    the solver with enforce_causality=False and central_pressure_pa explicitly.
    """

    epsilon_j_m3: float
    name: str = "Incompressible analytic validation (acausal)"
    validation_only: bool = True

    def __post_init__(self) -> None:
        if not isfinite(self.epsilon_j_m3) or self.epsilon_j_m3 <= 0:
            raise ValueError("epsilon must be positive")

    def h_from_rest_density(self, rho0: float) -> float:
        raise EOSDomainError("Incompressible test needs an independent central pressure")

    def h_from_pressure(self, pressure_pa: float) -> float:
        return log1p(_nonnegative(pressure_pa, "pressure")/self.epsilon_j_m3)

    def state_from_h(self, h: float) -> EOSState:
        return EOSState(self.epsilon_j_m3*expm1(_nonnegative(h,"enthalpy")),
                        self.epsilon_j_m3, float("nan"), float("inf"))

    def metadata(self) -> dict:
        return {"type":"incompressible", "name":self.name, "epsilon_j_m3":self.epsilon_j_m3,
                "status":"analytic validation ONLY, infinite sound speed"}


class TabulatedEOS(EquationOfState):
    """Strict cold, continuous 1D EOS import with no high-density extrapolation.

    A matched polytropic tail is required explicitly. At the first table row,
    epsilon=(1+a)rho0 c²+P/(Gamma_tail-1) fixes a; h=ln[(epsilon+P)/(rho0 c²(1+a))].
    This matches P, epsilon and chemical potential, not a physical crust model.
    Table thermodynamic errors and interpolation errors remain distinct from
    solve_ivp tolerance; a report is available as interpolation_diagnostics.
    """

    def __init__(self, pressure_pa: np.ndarray, energy_density_j_m3: np.ndarray,
                 rest_density_kg_m3: np.ndarray, *, tail_gamma: float,
                 provenance: str, enthalpy_grid_points: int = 8193,
                 first_law_tolerance: float = 0.05,
                 name: str = "Imported cold EOS", source_sha256: str | None = None) -> None:
        self.name, self.provenance, self.source_sha256 = name, provenance, source_sha256
        if not provenance.strip():
            raise ValueError("A source/provenance statement is required")
        arrays = [np.asarray(a, dtype=float).copy() for a in
                  (pressure_pa, energy_density_j_m3, rest_density_kg_m3)]
        if any(a.ndim != 1 for a in arrays) or len({len(a) for a in arrays}) != 1 or len(arrays[0]) < 8:
            raise ValueError("Need at least eight equal-length 1D rows")
        for a in arrays:
            if not np.all(np.isfinite(a)) or np.any(a <= 0) or np.any(np.diff(a) <= 0):
                raise EOSDomainError("P, epsilon, rho0 must be positive, finite and strictly increasing; density jumps unsupported")
        if (not isfinite(tail_gamma) or tail_gamma <= 1 or
            not 0 < first_law_tolerance < 1 or enthalpy_grid_points < 257):
            raise ValueError("Invalid tail, interpolation grid or thermodynamic tolerance")
        self.pressures, self.energies, self.densities = arrays
        self.tail_gamma, self.first_law_tolerance = tail_gamma, first_law_tolerance
        p0, e0, r0 = (a[0] for a in arrays)
        self.tail_K = p0/r0**tail_gamma
        self.tail_A = (e0-p0/(tail_gamma-1))/(r0*C*C)
        if self.tail_A <= 0:
            raise EOSDomainError("The chosen matched tail has a nonpositive rest-energy coefficient")
        self.h0 = log1p(tail_gamma/(tail_gamma-1)*p0/(self.tail_A*r0*C*C))
        self._loge = PchipInterpolator(np.log(self.pressures), np.log(self.energies), extrapolate=False)
        self._logrho = PchipInterpolator(np.log(self.pressures), np.log(self.densities), extrapolate=False)
        # Check the ACTUAL e(P),rho(P) interpolants used by the solver:
        # (d epsilon/d log P)/(d log rho/d log P) = epsilon+P.
        # Endpoints are included, with one-sided PCHIP derivatives.
        probe_logp = np.linspace(np.log(self.pressures[0]),np.log(self.pressures[-1]),max(1025,len(self.pressures)*5))
        e = np.exp(self._loge(probe_logp)); p = np.exp(probe_logp)
        de_dlp = self._loge.derivative()(probe_logp)*e
        dlrho_dlp = self._logrho.derivative()(probe_logp)
        if np.any(dlrho_dlp <= 0):
            raise EOSDomainError("Nonpositive density derivative in the interpolated table")
        residual = np.abs(de_dlp/(dlrho_dlp*(e+p))-1)
        self.first_law_max_relative_error = float(np.max(residual))
        if self.first_law_max_relative_error > first_law_tolerance:
            raise EOSDomainError(f"Table violates the cold first law: max residual {self.first_law_max_relative_error:.3g} > {first_law_tolerance}")
        # A uniform dense log-P grid avoids nearly coincident source/grid knots
        # that amplify roundoff in the derivative of the inverse PCHIP.
        # Source- and enthalpy-grid refinements are validated independently.
        grid = np.linspace(np.log(p0), np.log(self.pressures[-1]), enthalpy_grid_points)
        pp = np.exp(grid); ee = np.exp(self._loge(grid))
        hh = self.h0+cumulative_simpson(pp/(ee+pp), x=grid, initial=0.)
        self.max_enthalpy = float(hh[-1])
        self._h_of_logp = PchipInterpolator(grid, hh, extrapolate=False)
        self._logp_of_h = PchipInterpolator(hh, grid, extrapolate=False)
        probe = (hh[1:]+hh[:-1])/2
        lp = self._logp_of_h(probe)
        p = np.exp(lp); eps = np.exp(self._loge(lp))
        derivative = p*self._logp_of_h.derivative()(probe)
        self.enthalpy_identity_max_relative_error = float(np.max(np.abs(derivative/(eps+p)-1)))
        self.enthalpy_grid_points = int(len(grid))
        self.interpolation_diagnostics = {
            "cold_first_law_max_relative_error":self.first_law_max_relative_error,
            "enthalpy_identity_max_relative_error":self.enthalpy_identity_max_relative_error,
            "enthalpy_grid_points":self.enthalpy_grid_points,
            "warning":"These are interpolation diagnostics, not physical uncertainties.",
        }

    def _tail_state(self, h: float) -> EOSState:
        g, A = self.tail_gamma, self.tail_A
        rho = (A*expm1(h)*(g-1)*C*C/(g*self.tail_K))**(1/(g-1))
        p = self.tail_K*rho**g
        eps = A*rho*C*C+p/(g-1)
        q = g*p/(eps+p) if rho else 0.
        return EOSState(p, eps, rho, q)

    def state_from_h(self, h: float) -> EOSState:
        h = _nonnegative(h, "enthalpy")
        if h < self.h0:
            return self._tail_state(h)
        if h > self.max_enthalpy*(1+1e-13):
            raise EOSDomainError("Enthalpy exceeds the table; extrapolation is disabled")
        lp = float(self._logp_of_h(min(h, self.max_enthalpy)))
        p, eps = exp(lp), exp(float(self._loge(lp)))
        rho = exp(float(self._logrho(lp)))
        dloge_dlogp = float(self._loge.derivative()(lp))
        if dloge_dlogp <= 0:
            raise EOSDomainError("Nonpositive interpolated compressibility")
        return EOSState(p, eps, rho, p/(eps*dloge_dlogp))

    def h_from_rest_density(self, rho0: float) -> float:
        rho = _nonnegative(rho0, "rest density")
        if rho < self.densities[0]:
            g = self.tail_gamma
            return log1p(g*self.tail_K*rho**(g-1)/((g-1)*self.tail_A*C*C))
        if rho > self.densities[-1]*(1+1e-13):
            raise EOSDomainError("Density exceeds the table; extrapolation is disabled")
        target = log(min(rho, self.densities[-1]))
        # Invert the same interpolant, rather than interpolating a second,
        # approximately inverse map. This preserves the requested rho_c.
        lp = brentq(lambda v:float(self._logrho(v))-target,
                    log(self.pressures[0]),log(self.pressures[-1]),xtol=1e-13,rtol=1e-14)
        return float(self._h_of_logp(lp))

    def h_from_pressure(self, pressure_pa: float) -> float:
        p = _nonnegative(pressure_pa, "pressure")
        if p < self.pressures[0]:
            return self.h_from_rest_density((p/self.tail_K)**(1/self.tail_gamma))
        if p > self.pressures[-1]*(1+1e-13):
            raise EOSDomainError("Pressure exceeds the table; extrapolation is disabled")
        return float(self._h_of_logp(log(min(p, self.pressures[-1]))))

    @classmethod
    def from_csv(cls, path: str | Path, *, units: str, tail_gamma: float,
                 provenance: str, baryon_mass_kg: float | None = None, **kwargs) -> "TabulatedEOS":
        """Import a canonical CSV, not raw multidimensional CompOSE files.

        si: rho0_kg_m3,pressure_pa,epsilon_j_m3
        nuclear: n_b_fm3,pressure_mev_fm3,epsilon_mev_fm3 (m_b mandatory)
        cgs: rho0_g_cm3,pressure_dyn_cm2,epsilon_erg_cm3
        """
        path = Path(path)
        raw = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
        try:
            if units == "si":
                rho, p, eps = raw["rho0_kg_m3"], raw["pressure_pa"], raw["epsilon_j_m3"]
            elif units == "nuclear":
                if baryon_mass_kg is None or not isfinite(baryon_mass_kg) or baryon_mass_kg <= 0:
                    raise ValueError("A positive baryon_mass_kg convention is mandatory for nuclear units")
                rho = raw["n_b_fm3"]*1e45*baryon_mass_kg
                p = raw["pressure_mev_fm3"]*MEV_FM3_TO_PA
                eps = raw["epsilon_mev_fm3"]*MEV_FM3_TO_PA
            elif units == "cgs":
                rho = raw["rho0_g_cm3"]*1e3
                p, eps = raw["pressure_dyn_cm2"]*.1, raw["epsilon_erg_cm3"]*.1
            else:
                raise ValueError("units must be 'si', 'nuclear' or 'cgs'")
        except (ValueError, KeyError) as exc:
            raise ValueError(f"Invalid CSV schema or unit convention: {exc}") from exc
        table = cls(p, eps, rho, tail_gamma=tail_gamma, provenance=provenance,
                    source_sha256=sha256(path.read_bytes()).hexdigest(), **kwargs)
        table.input_units = units
        table.baryon_mass_kg = baryon_mass_kg
        return table

    def metadata(self) -> dict:
        return {"type":"table", "name":self.name, "provenance":self.provenance,
                "source_sha256":self.source_sha256, "rows":len(self.pressures),
                "tail_gamma":self.tail_gamma, "tail_rest_energy_coefficient":self.tail_A,
                "input_units":getattr(self,"input_units","si arrays"),
                "baryon_mass_kg":getattr(self,"baryon_mass_kg",None),
                "first_law_tolerance":self.first_law_tolerance,
                "interpolation":self.interpolation_diagnostics,
                "status":"source-dependent physical content; explicit phenomenological surface tail"}


def eos_from_spec(spec: dict[str, Any]) -> EquationOfState:
    """Serializable factory for CLI, saved experiments and process workers."""
    kind = spec["type"]
    if kind == "fermi":
        return DegenerateNeutronGas(neutron_mass_kg=spec.get("neutron_mass_kg", M_NEUTRON))
    if kind == "polytrope":
        return RelativisticPolytrope(K=spec["K_SI"], gamma=spec["gamma"])
    if kind == "piecewise":
        return PiecewisePolytrope(spec["K0_SI"], spec["gammas"], spec["breaks_kg_m3"])
    if kind == "css":
        return ConstantSoundSpeed(spec["s"], spec["epsilon_surface_j_m3"], spec["rho0_surface_kg_m3"])
    if kind == "table":
        return TabulatedEOS.from_csv(spec["path"], units=spec["units"],
            tail_gamma=spec["tail_gamma"], provenance=spec["provenance"],
            baryon_mass_kg=spec.get("baryon_mass_kg"))
    raise ValueError(f"Unknown EOS type: {kind}")
