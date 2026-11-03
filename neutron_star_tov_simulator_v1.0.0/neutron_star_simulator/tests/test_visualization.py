import tempfile
from pathlib import Path
import unittest
import numpy as np
from neutron_star_simulator.equations_of_state import DegenerateNeutronGas
from neutron_star_simulator.stability import sweep_central_densities
from neutron_star_simulator.tov_solver import solve_star
from neutron_star_simulator.visualization import plot_profiles,plot_sequences,plot_equations_of_state


class VisualizationTests(unittest.TestCase):
    def test_real_outputs_png_and_pdf(self):
        eos=DegenerateNeutronGas()
        with tempfile.TemporaryDirectory() as tmp:
            files=plot_profiles([solve_star(eos,1e18)],tmp)
            seq=sweep_central_densities(eos,np.geomspace(1e17,1e19,6))
            files+=plot_sequences({"Fermi":seq},tmp)
            files+=plot_equations_of_state({"Fermi":eos},tmp)
            self.assertEqual(len(files),18)
            for f in files:
                b=Path(f).read_bytes()
                self.assertGreater(len(b),5000)
                self.assertTrue(b.startswith(b"%PDF") if f.endswith(".pdf") else b.startswith(b"\x89PNG"))
