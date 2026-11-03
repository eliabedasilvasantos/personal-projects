from dataclasses import replace
import unittest
import numpy as np
from neutron_star_simulator.constants import C, G, GeometricUnits
from neutron_star_simulator.config import SolverConfig
from neutron_star_simulator.equations_of_state import (
    DegenerateNeutronGas, RelativisticPolytrope, IncompressibleEOS, ConstantSoundSpeed,
)
from neutron_star_simulator.tov_solver import integrate_tov, tov_rhs_si, geometric_rhs
from neutron_star_simulator.exceptions import CausalityError, SurfaceNotFoundError, IntegrationError


class SolverTests(unittest.TestCase):
    def test_si_and_geometric_equivalence(self):
        state = DegenerateNeutronGas().state_from_rest_density(1e18)
        r, m, h = 8500., 1e30, .1
        u = GeometricUnits(5000)
        rhs = geometric_rhs(r/u.length_m,np.array([m/u.mass_kg,h,0.]),state,u)
        dm, dp = tov_rhs_si(r,m,state.pressure_pa,state.energy_density_j_m3)
        self.assertAlmostEqual(rhs[0]*u.mass_kg/u.length_m/dm,1,places=14)
        self.assertAlmostEqual(rhs[1]/u.length_m*(state.energy_density_j_m3+state.pressure_pa)/dp,1,places=14)

    def test_incompressible_exact(self):
        rho, R = 7e17, 11_000.
        eps = rho*C*C; M = 4*np.pi/3*rho*R**3; compactness = G*M/(R*C*C)
        a = np.sqrt(1-2*compactness)
        pc = eps*(1-a)/(3*a-1)
        cfg = SolverConfig(rtol=1e-10,atol_mass=1e-13,atol_enthalpy=1e-13,
                           atol_baryonic_mass=1e-13,enforce_causality=False)
        s = integrate_tov(IncompressibleEOS(eps),config=cfg,central_pressure_pa=pc)
        self.assertLess(abs(s.radius_m[-1]/R-1),1e-8)
        self.assertLess(abs(s.gravitational_mass_kg[-1]/M-1),1e-8)
        b = np.sqrt(1-2*compactness*(s.radius_m/R)**2)
        exact_p = eps*(b-a)/(3*a-b)
        self.assertLess(np.max(np.abs(s.pressure_pa-exact_p))/pc,1e-8)

    def test_center_and_surface(self):
        s = integrate_tov(DegenerateNeutronGas(),1e18)
        self.assertEqual(s.radius_m[0],0)
        self.assertEqual(s.gravitational_mass_kg[0],0)
        self.assertEqual(s.pressure_pa[-1],0)
        self.assertLess(abs(s.diagnostics["surface_enthalpy_residual"]),1e-12)
        self.assertTrue(np.all(np.diff(s.gravitational_mass_kg)>=-1e-9*s.gravitational_mass_kg[-1]))
        self.assertTrue(np.all(np.diff(s.pressure_pa)<=0))
        self.assertTrue(np.all(np.isfinite(s.energy_density_j_m3)))

    def test_causality_rejects(self):
        acausal = RelativisticPolytrope.from_reference(1e18,1e36,3.)
        with self.assertRaises(CausalityError): integrate_tov(acausal,1e18)
        with self.assertRaises(CausalityError):
            integrate_tov(IncompressibleEOS(1e35),central_pressure_pa=1e34)

    def test_no_surface_is_failure(self):
        with self.assertRaises(SurfaceNotFoundError):
            integrate_tov(DegenerateNeutronGas(),1e18,SolverConfig(max_radius_m=1000))

    def test_large_central_start_rejects(self):
        with self.assertRaises(IntegrationError):
            integrate_tov(DegenerateNeutronGas(),1e18,SolverConfig(initial_radius_m=5000))

    def test_no_ambiguous_central_condition(self):
        for kwargs in ({}, {"central_rest_density_kg_m3":1e18,"central_pressure_pa":1e34}):
            with self.assertRaises(ValueError): integrate_tov(DegenerateNeutronGas(),**kwargs)

    def test_method_and_length_scale_invariance(self):
        base = integrate_tov(DegenerateNeutronGas(),1e18)
        for cfg in (SolverConfig(method="RK45",rtol=1e-10,atol_enthalpy=1e-13,atol_mass=1e-13),
                    SolverConfig(length_scale_m=1000), SolverConfig(length_scale_m=100_000)):
            other = integrate_tov(DegenerateNeutronGas(),1e18,cfg)
            self.assertLess(abs(other.radius_m[-1]/base.radius_m[-1]-1),2e-6)
            self.assertLess(abs(other.gravitational_mass_kg[-1]/base.gravitational_mass_kg[-1]-1),2e-6)

    def test_self_bound_surface(self):
        eos = ConstantSoundSpeed(1/3,4e34,4e34/C**2)
        s = integrate_tov(eos,1e18)
        self.assertEqual(s.pressure_pa[-1],0)
        self.assertEqual(s.energy_density_j_m3[-1],4e34)
        self.assertGreater(s.radius_m[-1],1000)
