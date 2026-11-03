import json
import tempfile
from pathlib import Path
import unittest
import numpy as np
from scipy.integrate import simpson
from neutron_star_simulator.constants import C, G, M_SUN_KG
from neutron_star_simulator.equations_of_state import DegenerateNeutronGas, IncompressibleEOS
from neutron_star_simulator.config import SolverConfig
from neutron_star_simulator.tov_solver import solve_star


class PropertyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.star = solve_star(DegenerateNeutronGas(),1e18,SolverConfig(profile_points=2001))

    def test_mass_integral(self):
        s = self.star.solution
        m = simpson(4*np.pi*s.radius_m**2*s.energy_density_j_m3/C**2,x=s.radius_m)
        self.assertLess(abs(m/self.star.mass_kg-1),2e-6)

    def test_baryonic_mass_proper_volume(self):
        s = self.star.solution
        proper = np.sqrt(1-2*self.star.compactness_profile)
        mb = simpson(4*np.pi*s.radius_m**2*s.rest_density_kg_m3/proper,x=s.radius_m)
        self.assertLess(abs(mb/self.star.baryonic_mass_kg-1),2e-6)
        self.assertGreater(self.star.binding_energy_j,0)

    def test_redshift_and_lapse(self):
        self.assertAlmostEqual((1+self.star.surface_redshift)*self.star.lapse_profile[-1],1,places=14)
        self.assertEqual(self.star.compactness_profile[0],0)
        self.assertAlmostEqual(self.star.compactness,G*self.star.mass_kg/(self.star.radius_m*C*C))
        self.assertTrue(np.all(np.diff(self.star.lapse_profile)>=0))
        self.assertLess(np.ptp(np.log(self.star.lapse_profile)+self.star.solution.enthalpy),1e-14)

    def test_baryonic_mass_unavailable_for_incompressible(self):
        star = solve_star(IncompressibleEOS(7e17*C*C),config=SolverConfig(enforce_causality=False),
                          central_pressure_pa=1e34)
        self.assertIsNone(star.baryonic_mass_kg)
        self.assertIsNone(star.binding_energy_j)
        json.dumps(star.summary(),allow_nan=False)

    def test_export_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = self.star.export(tmp)
            data = np.genfromtxt(files["profile_csv"],delimiter=",",names=True)
            summary = json.loads(Path(files["summary_json"]).read_text())
            self.assertEqual(len(data),2001)
            self.assertAlmostEqual(data["gravitational_mass_kg"][-1]/M_SUN_KG,self.star.mass_solar)
            self.assertEqual(summary["radius_m"],self.star.radius_m)
            self.assertNotIn("NaN",Path(files["summary_json"]).read_text())
