"""Audited SI constants and explicit conversions (CODATA 2022; IAU 2015 B3).

The solar *mass parameter* is nominal/exact. M_SUN_KG = GM_SUN/G is a
conversion convention, not an exact measurement of the Sun's mass. G and
the neutron mass have measurement uncertainty; uncertainties are not propagated.
"""

from dataclasses import dataclass
from math import pi

G = 6.67430e-11  # m^3 kg^-1 s^-2, CODATA 2022
C = 299_792_458.0  # m s^-1, exact
PLANCK = 6.62607015e-34  # J s, exact
HBAR = PLANCK / (2.0 * pi)
ELEMENTARY_CHARGE = 1.602176634e-19  # C; eV -> J, exact
M_NEUTRON = 1.67492750056e-27  # kg, CODATA 2022
GM_SUN = 1.3271244e20  # m^3 s^-2, IAU nominal conversion constant
M_SUN_KG = GM_SUN / G
MEV_FM3_TO_PA = 1e6 * ELEMENTARY_CHARGE / 1e-45
FM3_INV_TO_M3_INV = 1e45
G_CM3_TO_KG_M3 = 1e3
DYN_CM2_TO_PA = 0.1


@dataclass(frozen=True)
class GeometricUnits:
    """Dimensionless geometrized variables based on one length scale L.

    x=r/L, mu=G*m/(c^2 L), e=G*L^2*epsilon/c^4, p=G*L^2*P/c^4.
    rho0_dimensionless=G*L^2*rho0/c^2. No implicit c=1 conversion of kg/m^3.
    """

    length_m: float = 10_000.0

    def __post_init__(self) -> None:
        from math import isfinite
        if not isfinite(self.length_m) or self.length_m <= 0:
            raise ValueError("length_m must be positive and finite")

    @property
    def mass_kg(self) -> float:
        return C**2 * self.length_m / G

    @property
    def pressure_pa(self) -> float:
        return C**4 / (G * self.length_m**2)

    @property
    def density_kg_m3(self) -> float:
        return C**2 / (G * self.length_m**2)


def mev_fm3_to_si(value: float) -> float:
    """Convert either pressure or total energy density from MeV/fm^3 to Pa=J/m^3."""
    return value * MEV_FM3_TO_PA


def baryon_density_to_rest_density(n_fm3: float, baryon_mass_kg: float) -> float:
    """Convert n_b in fm^-3 to rho0 in kg/m^3 with an EXPLICIT mass convention."""
    return n_fm3 * FM3_INV_TO_M3_INV * baryon_mass_kg
