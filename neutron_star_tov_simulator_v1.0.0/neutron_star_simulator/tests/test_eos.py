import tempfile
from pathlib import Path
import unittest
import numpy as np
from scipy.integrate import quad
from neutron_star_simulator.constants import C, M_NEUTRON, HBAR, MEV_FM3_TO_PA
from neutron_star_simulator.equations_of_state import (
    DegenerateNeutronGas, RelativisticPolytrope, PiecewisePolytrope,
    ConstantSoundSpeed, TabulatedEOS, eos_from_spec,
)
from neutron_star_simulator.exceptions import EOSDomainError


class EOSTests(unittest.TestCase):
    def setUp(self):
        self.poly = RelativisticPolytrope(.01, 2.)
        self.fermi = DegenerateNeutronGas()
        self.piece = PiecewisePolytrope(.01, [2, 2.3, 1.8], [2e17, 8e17])
        self.css = ConstantSoundSpeed(1/3, 4e34, 4e34/C**2)
        self.models = [self.poly, self.fermi, self.piece, self.css]

    def test_density_and_pressure_roundtrip(self):
        for eos in self.models:
            for rho in np.geomspace(5e17, 5e18, 7):
                with self.subTest(eos=eos.name, rho=rho):
                    h = eos.h_from_rest_density(rho)
                    state = eos.state_from_h(h)
                    self.assertAlmostEqual(state.rest_density_kg_m3/rho, 1, places=11)
                    self.assertAlmostEqual(eos.h_from_pressure(state.pressure_pa)/h, 1, places=10)

    def test_cold_first_law(self):
        for eos in self.models:
            for rho in (5e17, 1e18, 3e18):
                d = 1e-5
                lo, hi = (eos.state_from_rest_density(rho*np.exp(s*d)) for s in (-1,1))
                mid = eos.state_from_rest_density(rho)
                derivative = (hi.energy_density_j_m3-lo.energy_density_j_m3)/(2*d)
                self.assertLess(abs(derivative/(mid.energy_density_j_m3+mid.pressure_pa)-1), 2e-8)

    def test_enthalpy_and_sound_speed_derivatives(self):
        for eos in self.models:
            h = eos.h_from_rest_density(1e18); d = h*1e-5
            lo, hi, mid = (eos.state_from_h(v) for v in (h-d, h+d, h))
            dp = hi.pressure_pa-lo.pressure_pa
            de = hi.energy_density_j_m3-lo.energy_density_j_m3
            self.assertLess(abs(dp/(2*d)/(mid.energy_density_j_m3+mid.pressure_pa)-1), 2e-8)
            self.assertLess(abs(dp/de/mid.sound_speed_squared_over_c2-1), 2e-8)

    def test_fermi_independent_momentum_quadrature(self):
        B = M_NEUTRON**4*C**5/(np.pi**2*HBAR**3)
        for x in (1e-5, .01, .079999, .080001, .1, 1., 20.):
            s = self.fermi.state_from_h(.5*np.log1p(x*x))
            ie = x**3*quad(lambda u: u*u*np.sqrt(1+(x*u)**2), 0, 1, epsabs=1e-13)[0]
            ip = x**5/3*quad(lambda u: u**4/np.sqrt(1+(x*u)**2), 0, 1, epsabs=1e-13)[0]
            self.assertLess(abs(s.energy_density_j_m3/(B*ie)-1), 2e-10)
            self.assertLess(abs(s.pressure_pa/(B*ip)-1), 2e-10)

    def test_fermi_limits(self):
        low = self.fermi.state_from_h(.5*np.log1p(1e-8))
        self.assertLess(abs(low.energy_density_j_m3/(low.rest_density_kg_m3*C*C)-1), 1e-7)
        high = self.fermi.state_from_h(.5*np.log1p(1e10))
        self.assertLess(abs(high.pressure_pa/high.energy_density_j_m3-1/3), 1e-9)

    def test_piecewise_continuity(self):
        for rho in self.piece.breaks:
            lo, hi = (self.piece.state_from_rest_density(rho*(1+s*1e-10)) for s in (-1,1))
            for key in ("pressure_pa","energy_density_j_m3","rest_density_kg_m3"):
                self.assertLess(abs(getattr(hi,key)/getattr(lo,key)-1), 1e-8)

    def test_invalid_and_surface_states(self):
        for eos in self.models:
            with self.assertRaises(EOSDomainError): eos.state_from_h(-.1)
        self.assertEqual(self.poly.state_from_h(0).pressure_pa, 0)
        self.assertEqual(self.fermi.state_from_h(0).energy_density_j_m3, 0)
        self.assertEqual(self.css.state_from_h(0).energy_density_j_m3, 4e34)
        with self.assertRaises(EOSDomainError): self.css.h_from_rest_density(1e10)
        with self.assertRaises(ValueError): RelativisticPolytrope(.01, 1)
        with self.assertRaises(ValueError): ConstantSoundSpeed(2, 1e34, 1e17)

    def make_table(self, **kwargs):
        rho = np.geomspace(1e10, 1e19, 301)
        return TabulatedEOS(.01*rho**2, rho*C*C+.01*rho**2, rho,
                            tail_gamma=2., provenance="SYNTHETIC test polytrope", **kwargs)

    def test_table_roundtrip_and_diagnostics(self):
        table = self.make_table()
        for rho in np.geomspace(1e8, 8e18, 17):
            exact = self.poly.state_from_rest_density(rho)
            actual = table.state_from_rest_density(rho)
            self.assertLess(abs(actual.pressure_pa/exact.pressure_pa-1), 2e-7)
            self.assertLess(abs(actual.energy_density_j_m3/exact.energy_density_j_m3-1), 3e-6)
        self.assertLess(table.first_law_max_relative_error, 1e-3)
        self.assertLess(table.enthalpy_identity_max_relative_error, 1e-4)
        with self.assertRaises(EOSDomainError): table.h_from_rest_density(2e19)
        with self.assertRaises(EOSDomainError): table.state_from_h(table.max_enthalpy*1.1)

    def test_table_rejects_bad_units_or_first_law(self):
        rho = np.geomspace(1e10, 1e19, 101)
        with self.assertRaises(EOSDomainError):
            TabulatedEOS(.01*rho**2, rho*C*C*(rho/1e15)**.1+.01*rho**2, rho,
                         tail_gamma=2, provenance="bad", first_law_tolerance=.05)
        with self.assertRaises(EOSDomainError):
            TabulatedEOS(np.ones(101),rho*C*C,rho,tail_gamma=2,provenance="bad")
        with self.assertRaises(ValueError):
            TabulatedEOS(.01*rho**2,rho*C*C+.01*rho**2,rho,tail_gamma=2,provenance="")

    def test_csv_nuclear_conversion(self):
        rho = np.geomspace(1e10,1e19,101); p = .01*rho*rho; eps = rho*C*C+p
        with tempfile.TemporaryDirectory() as tmp:
            f = Path(tmp)/"eos.csv"
            np.savetxt(f,np.c_[rho/(M_NEUTRON*1e45), p/MEV_FM3_TO_PA, eps/MEV_FM3_TO_PA],
                       delimiter=",",header="n_b_fm3,pressure_mev_fm3,epsilon_mev_fm3",comments="")
            table = TabulatedEOS.from_csv(f,units="nuclear",tail_gamma=2,provenance="synthetic",baryon_mass_kg=M_NEUTRON)
            self.assertAlmostEqual(table.densities[-1]/rho[-1],1,places=12)
            self.assertEqual(len(table.source_sha256),64)
            with self.assertRaises(ValueError):
                TabulatedEOS.from_csv(f,units="nuclear",tail_gamma=2,provenance="synthetic")

    def test_factory_roundtrip(self):
        for eos in self.models:
            other = eos_from_spec(eos.metadata())
            self.assertAlmostEqual(other.state_from_rest_density(1e18).pressure_pa /
                                   eos.state_from_rest_density(1e18).pressure_pa,1,places=12)
