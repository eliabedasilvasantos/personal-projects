import unittest
from neutron_star_simulator.validation import (
    dimensional_analysis,incompressible_benchmark,newtonian_n1_benchmark,
    fermi_published_benchmark,convergence_study,table_interpolation_study,profile_consistency,
)
from neutron_star_simulator.equations_of_state import DegenerateNeutronGas
from neutron_star_simulator.tov_solver import solve_star


class ValidationTests(unittest.TestCase):
    def test_dimensions(self): self.assertTrue(dimensional_analysis()["passed"])
    def test_exact_solution(self): self.assertTrue(incompressible_benchmark()["passed"])
    def test_newtonian_asymptotic_limit(self): self.assertTrue(newtonian_n1_benchmark()["passed"])
    def test_published_table(self): self.assertTrue(fermi_published_benchmark()["passed"])
    def test_tolerance_center_scale_and_step(self): self.assertTrue(convergence_study()["passed"])
    def test_table_discretization(self): self.assertTrue(table_interpolation_study()["passed"])
    def test_independent_profile_checks(self):
        self.assertTrue(profile_consistency(solve_star(DegenerateNeutronGas(),1e18))["passed"])
