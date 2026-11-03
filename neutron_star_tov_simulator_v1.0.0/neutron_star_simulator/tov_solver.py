"""Adaptive TOV integration in areal radius and relativistic enthalpy.

The physical equations are exposed in SI for dimensional/equivalence tests.
The actual solve_ivp state is [mu, h, mu_b], dimensionless. No pressure floor
or Newtonian approximation is used in the integration.
"""

from dataclasses import dataclass
from math import isfinite, pi, sqrt
from time import perf_counter
from typing import Any
import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from .config import SolverConfig
from .constants import C, G, GeometricUnits
from .equations_of_state import EquationOfState, EOSState
from .exceptions import CausalityError, EOSDomainError, IntegrationError, SurfaceNotFoundError

FloatArray = NDArray[np.float64]


def tov_rhs_si(radius_m: float, mass_kg: float, pressure_pa: float,
               energy_density_j_m3: float) -> tuple[float, float]:
    """Exact TOV RHS in kg/m and Pa/m; center must use a regular series."""
    if radius_m <= 0 or mass_kg < 0 or pressure_pa < 0 or energy_density_j_m3 < 0:
        raise ValueError("TOV SI requires r>0 and nonnegative matter/mass")
    metric = 1-2*G*mass_kg/(radius_m*C*C)
    if metric <= 0:
        raise IntegrationError("Trapped geometry: 1-2Gm/(rc²)<=0")
    dm = 4*pi*radius_m**2*energy_density_j_m3/C**2
    dp = -G*(energy_density_j_m3+pressure_pa)*(mass_kg+4*pi*radius_m**3*pressure_pa/C**2)/(
        C**2*radius_m**2*metric)
    return dm, dp


def geometric_rhs(x: float, y: FloatArray, state: EOSState,
                  units: GeometricUnits, baryons_available: bool = True) -> FloatArray:
    """Dimensionless (dmu/dx, dh/dx, dmu_b/dx) at x>0."""
    mu = float(y[0])
    metric = 1-2*mu/x
    if not isfinite(metric) or metric <= 0:
        raise IntegrationError("Trapped or nonfinite geometry encountered during an integration stage")
    e, p = state.energy_density_j_m3/units.pressure_pa, state.pressure_pa/units.pressure_pa
    dm = 4*pi*x*x*e
    dh = -(mu+4*pi*x**3*p)/(x*x*metric)
    db = (4*pi*x*x*state.rest_density_kg_m3/units.density_kg_m3/sqrt(metric)
          if baryons_available else 0.0)
    return np.array([dm, dh, db], dtype=float)


@dataclass(frozen=True)
class TOVSolution:
    radius_m: FloatArray
    gravitational_mass_kg: FloatArray
    enthalpy: FloatArray
    pressure_pa: FloatArray
    energy_density_j_m3: FloatArray
    rest_density_kg_m3: FloatArray
    baryonic_mass_kg: FloatArray
    sound_speed_squared_over_c2: FloatArray
    central_enthalpy: float
    central_pressure_pa: float
    central_energy_density_j_m3: float
    central_rest_density_kg_m3: float
    config: SolverConfig
    eos_metadata: dict[str, Any]
    diagnostics: dict[str, Any]


def _check_state(state: EOSState, config: SolverConfig) -> None:
    if not all(isfinite(x) and x >= 0 for x in
               (state.pressure_pa, state.energy_density_j_m3)):
        raise EOSDomainError("EOS returned nonfinite or negative pressure/energy")
    q = state.sound_speed_squared_over_c2
    if q < 0 or (state.pressure_pa > 0 and q == 0) or np.isnan(q):
        raise CausalityError("Nonpositive compressibility or invalid sound speed")
    if config.enforce_causality and q > 1+1e-10:
        raise CausalityError(f"Acausal EOS state: c_s²/c²={q:.8g}; configuration rejected")


def integrate_tov(eos: EquationOfState,
                  central_rest_density_kg_m3: float | None = None,
                  config: SolverConfig | None = None, *,
                  central_pressure_pa: float | None = None) -> TOVSolution:
    """Compute one star, requiring exactly one central density or pressure.

    All physical outputs are SI. An equilibrium is returned only after the
    zero-enthalpy event, successful integration and profile sanity checks.
    Acausal analytic validation requires enforce_causality=False explicitly.
    """
    cfg = config or SolverConfig()
    if (central_rest_density_kg_m3 is None) == (central_pressure_pa is None):
        raise ValueError("Provide exactly one central rest density or central pressure")
    hc = (eos.h_from_rest_density(central_rest_density_kg_m3)
          if central_pressure_pa is None else eos.h_from_pressure(central_pressure_pa))
    if not isfinite(hc) or hc <= 0:
        raise EOSDomainError("The central enthalpy must be positive and finite")
    central = eos.state_from_h(hc)
    _check_state(central, cfg)
    if central.energy_density_j_m3 <= 0 or central.pressure_pa <= 0:
        raise EOSDomainError("A nonvacuum center with positive pressure is required")
    baryons = isfinite(central.rest_density_kg_m3)
    if baryons and central.rest_density_kg_m3 <= 0:
        raise EOSDomainError("Invalid central baryonic rest density")
    units = GeometricUnits(cfg.length_scale_m)
    x0 = cfg.initial_radius_m/units.length_m
    ec = central.energy_density_j_m3/units.pressure_pa
    pc = central.pressure_pa/units.pressure_pa
    bc = central.rest_density_kg_m3/units.density_kg_m3 if baryons else 0.0
    A = 2*pi*(ec/3+pc)
    correction = A*x0*x0
    if correction/hc > 1e-3:
        raise IntegrationError("Initial radius too large for the central series; reduce initial_radius_m")
    y0 = np.array([4*pi/3*ec*x0**3, hc-correction, 4*pi/3*bc*x0**3])
    max_sound = central.sound_speed_squared_over_c2

    def rhs(x: float, y: FloatArray) -> FloatArray:
        nonlocal max_sound
        # Only the interior contributes to output. Trial steps may cross h=0;
        # extending the EOS at its surface prevents undefined fractional powers.
        state = eos.state_from_h(max(float(y[1]), 0.0))
        _check_state(state, cfg)
        max_sound = max(max_sound, state.sound_speed_squared_over_c2)
        return geometric_rhs(x, y, state, units, baryons)

    def surface(x: float, y: FloatArray) -> float:
        return float(y[1])
    surface.terminal = True
    surface.direction = -1.0

    def horizon(x: float, y: FloatArray) -> float:
        return 1-2*float(y[0])/x-cfg.horizon_margin
    horizon.terminal = True
    horizon.direction = -1.0

    started = perf_counter()
    sol = solve_ivp(rhs, (x0, cfg.max_radius_m/units.length_m), y0,
                    method=cfg.method, rtol=cfg.rtol,
                    atol=[cfg.atol_mass, cfg.atol_enthalpy, cfg.atol_baryonic_mass],
                    max_step=cfg.max_step_m/units.length_m,
                    events=(surface, horizon), dense_output=True)
    if not sol.success:
        raise IntegrationError(f"solve_ivp failed: {sol.message}")
    if len(sol.t_events[1]):
        raise IntegrationError("Near-horizon termination; no stellar surface returned")
    if not len(sol.t_events[0]):
        raise SurfaceNotFoundError(f"P=0 not reached by {cfg.max_radius_m/1000:g} km; enlarge max_radius_m or inspect the EOS")
    xs = float(sol.t_events[0][0])
    raw_surface_h = float(sol.y_events[0][0,1])
    x = np.r_[0., np.linspace(x0, xs, cfg.profile_points-1)]
    y = sol.sol(x[1:])
    h = np.r_[hc, y[1]]
    # Only set the root's roundoff residual to its known boundary value.
    if np.min(h[:-1]) < 0:
        raise IntegrationError("Negative interior enthalpy in dense output")
    h[-1] = 0.0
    states = [eos.state_from_h(float(v)) for v in h]
    for state in states:
        _check_state(state, cfg)
    m = np.r_[0., y[0]]*units.mass_kg
    mb = np.r_[0., y[2]]*units.mass_kg if baryons else np.full_like(m, np.nan)
    p = np.array([s.pressure_pa for s in states])
    e = np.array([s.energy_density_j_m3 for s in states])
    rho = np.array([s.rest_density_kg_m3 for s in states])
    q = np.array([s.sound_speed_squared_over_c2 for s in states])
    if np.any(np.diff(m) < -1e-9*m[-1]) or np.any(np.diff(h) > 1e-9*hc):
        raise IntegrationError("Nonmonotonic mass or enthalpy profile; tighten integration settings")
    r = x*units.length_m
    if m[-1] <= 0 or r[-1] <= 0 or 2*G*m[-1]/(r[-1]*C*C) >= 1:
        raise IntegrationError("Invalid total mass, radius or surface geometry")
    for array in (r, m, h, mb, p, e, rho, q):
        array.setflags(write=False)
    return TOVSolution(r, m, h, p, e, rho, mb, q, hc,
        central.pressure_pa, central.energy_density_j_m3, central.rest_density_kg_m3,
        cfg, eos.metadata(), {"success":True, "termination":"zero_enthalpy_surface",
        "surface_enthalpy_residual":raw_surface_h, "function_evaluations":sol.nfev,
        "accepted_steps":len(sol.t)-1, "message":sol.message,
        "wall_seconds":perf_counter()-started,
        "max_sound_speed_squared_over_c2_evaluated":max_sound,
        "central_series_relative_h_correction":correction/hc,
        "baryonic_mass_available":baryons,
        "causality_enforced":cfg.enforce_causality,
        "scientific_status":("analytic_validation_acausal" if eos.validation_only else
                             "hydrostatic_model_no_dynamic_mode_calculation")})


def solve_star(eos: EquationOfState, central_rest_density_kg_m3: float | None = None,
               config: SolverConfig | None = None, *, central_pressure_pa: float | None = None):
    """Convenient high-level wrapper with observables and CSV/JSON export."""
    from .stellar_properties import StellarModel
    return StellarModel(integrate_tov(eos, central_rest_density_kg_m3, config,
                                      central_pressure_pa=central_pressure_pa))
