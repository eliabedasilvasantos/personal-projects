import unittest
from neutron_star_simulator.constants import (
    C, G, GM_SUN, M_SUN_KG, MEV_FM3_TO_PA, GeometricUnits,
    baryon_density_to_rest_density, M_NEUTRON,
)
from neutron_star_simulator.config import SolverConfig


class UnitTests(unittest.TestCase):
    def test_energy_conversion(self):
        self.assertAlmostEqual(MEV_FM3_TO_PA / 1.602176634e32, 1.0, places=15)

    def test_solar_convention(self):
        self.assertAlmostEqual(G * M_SUN_KG / GM_SUN, 1.0, places=15)

    def test_geometric_roundtrip(self):
        for length in (1e3, 1e4, 1e5):
            u = GeometricUnits(length)
            mass, pressure, density = 2.7e30, 2.1e34, 8e17
            self.assertAlmostEqual((mass / u.mass_kg) * u.mass_kg / mass, 1., places=15)
            self.assertAlmostEqual((pressure / u.pressure_pa) * u.pressure_pa / pressure, 1., places=15)
            self.assertAlmostEqual(u.pressure_pa / (C*C * u.density_kg_m3), 1., places=15)
            self.assertAlmostEqual(G * density * length**2 / C**2, density / u.density_kg_m3)

    def test_baryon_density(self):
        self.assertAlmostEqual(baryon_density_to_rest_density(.16, M_NEUTRON) /
                               (.16e45*M_NEUTRON), 1., places=15)

    def test_invalid_settings(self):
        for kwargs in ({"rtol":-1}, {"initial_radius_m":1e9}, {"profile_points":3},
                       {"method":"unknown"}, {"length_scale_m":float("nan")}):
            with self.assertRaises(ValueError):
                SolverConfig(**kwargs)
