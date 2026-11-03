"""Numerical settings. Absolute tolerances apply to dimensionless state variables."""

from dataclasses import dataclass, asdict
from math import isfinite


@dataclass(frozen=True)
class SolverConfig:
    method: str = "DOP853"
    rtol: float = 1e-8
    atol_mass: float = 1e-11
    atol_enthalpy: float = 1e-11
    atol_baryonic_mass: float = 1e-11
    initial_radius_m: float = 0.5
    max_radius_m: float = 100_000.0
    max_step_m: float = 150.0
    length_scale_m: float = 10_000.0
    profile_points: int = 701
    enforce_causality: bool = True
    horizon_margin: float = 1e-6

    def __post_init__(self) -> None:
        if self.method not in {"DOP853", "RK45", "Radau", "BDF", "LSODA"}:
            raise ValueError("Unsupported solve_ivp method")
        for key in ("rtol", "atol_mass", "atol_enthalpy", "atol_baryonic_mass",
                    "initial_radius_m", "max_radius_m", "max_step_m", "length_scale_m"):
            if not isfinite(getattr(self, key)) or getattr(self, key) <= 0:
                raise ValueError(f"{key} must be positive and finite")
        if not 1e-13 <= self.rtol < 1:
            raise ValueError("rtol must be in [1e-13, 1)")
        if self.initial_radius_m >= self.max_radius_m:
            raise ValueError("The initial radius must be below the outer boundary")
        if not isinstance(self.profile_points, int) or self.profile_points < 51:
            raise ValueError("profile_points must be an integer >=51")
        if not 0 < self.horizon_margin < 0.1:
            raise ValueError("horizon_margin must be in (0, 0.1)")

    def to_dict(self) -> dict:
        return asdict(self)
