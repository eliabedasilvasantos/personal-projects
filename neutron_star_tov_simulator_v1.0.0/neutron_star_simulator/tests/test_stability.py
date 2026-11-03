import unittest
import numpy as np
from neutron_star_simulator.config import SolverConfig
from neutron_star_simulator.equations_of_state import DegenerateNeutronGas, RelativisticPolytrope
from neutron_star_simulator.stability import sweep_central_densities, analyze_stability, refine_maximum


class StabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.eos=DegenerateNeutronGas()
        cls.sequence=sweep_central_densities(cls.eos,np.geomspace(1e16,2e19,35))

    def test_physical_maximum_bracket_and_refine(self):
        analysis=analyze_stability(self.sequence,assume_low_density_stable=True)
        self.assertTrue(analysis["largest_sample_is_bracketed_peak"])
        self.assertIn("initial_branch_candidate_stable",analysis["labels"])
        self.assertIn("negative_slope_candidate_unstable",analysis["labels"])
        max_result=refine_maximum(self.eos,self.sequence)
        self.assertEqual(max_result["status"],"refined_local_maximum")
        m=max_result["model"]["gravitational_mass_solar"]
        self.assertTrue(.70<m<.72)

    def test_no_unstated_stability_assumption(self):
        labels=analyze_stability(self.sequence)["labels"]
        self.assertNotIn("initial_branch_candidate_stable",labels)
        self.assertIn("positive_slope_unclassified",labels)

    def test_boundary_not_maximum(self):
        seq=sweep_central_densities(self.eos,np.geomspace(1e16,1e17,5))
        result=refine_maximum(self.eos,seq)
        self.assertEqual(result["status"],"boundary_limited_no_maximum_claim")
        self.assertFalse(analyze_stability(seq)["largest_sample_is_bracketed_peak"])

    def test_failed_points_remain(self):
        eos=RelativisticPolytrope.from_reference(1e18,2e34,3)
        seq=sweep_central_densities(eos,np.geomspace(1e17,1e20,13))
        self.assertEqual(len(seq.points),13)
        self.assertFalse(seq.is_complete)
        a=analyze_stability(seq)
        for i,p in enumerate(seq.points):
            if p.model is None:
                self.assertEqual(a["labels"][i],"failed")
                self.assertIsNone(a["dM_depsilon_c_kg_m3_per_j"][i])
                self.assertIsNotNone(p.error_type)

    def test_rejects_unsorted_duplicates(self):
        for grid in ([1e18,1e17],[1e18,1e18],[-1,2],[]):
            with self.assertRaises(ValueError): sweep_central_densities(self.eos,grid)

    def test_density_grid_refinement_peak(self):
        coarse=refine_maximum(self.eos,self.sequence)
        fine=refine_maximum(self.eos,sweep_central_densities(self.eos,np.geomspace(1e16,2e19,60)))
        self.assertLess(abs(fine["model"]["gravitational_mass_solar"]/coarse["model"]["gravitational_mass_solar"]-1),1e-6)
