import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile
import numpy as np
from neutron_star_simulator.equations_of_state import DegenerateNeutronGas
from neutron_star_simulator.stability import sweep_central_densities
from neutron_star_simulator.tov_solver import solve_star
from neutron_star_simulator.streamlit_app import result_zip


class InterfaceTests(unittest.TestCase):
    def test_cli_single_equilibrium(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=subprocess.run([sys.executable,"-m","neutron_star_simulator.main","star",
                "--rho-c","1e18","--no-plots","--out",tmp],capture_output=True,text=True,check=False)
            self.assertEqual(result.returncode,0,result.stderr)
            report=json.loads((Path(tmp)/"star_summary.json").read_text())
            manifest=json.loads((Path(tmp)/"manifest.json").read_text())
            self.assertTrue(.62<report["gravitational_mass_solar"]<.64)
            self.assertIn("tov_solver.py",manifest["source_sha256"])
            self.assertEqual(len(manifest["source_sha256"]["tov_solver.py"]),64)

    def test_cli_rejects_unstated_polytrope_parameters(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=subprocess.run([sys.executable,"-m","neutron_star_simulator.main","star",
                "--eos","polytrope","--rho-c","1e18","--no-plots","--out",tmp],
                capture_output=True,text=True,check=False)
            self.assertEqual(result.returncode,1)
            self.assertIn("Explicit EOS parameters required",result.stderr)
            self.assertFalse((Path(tmp)/"star_summary.json").exists())

    def test_cli_failed_sequence_returns_nonzero_with_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=subprocess.run([sys.executable,"-m","neutron_star_simulator.main","sweep",
                "--rho-min","1e17","--rho-max","1e19","--points","8","--max-radius-km","1",
                "--no-plots","--out",tmp],capture_output=True,text=True,check=False)
            self.assertEqual(result.returncode,2,result.stderr)
            report=json.loads((Path(tmp)/"sequence.json").read_text())
            self.assertEqual(len(report["points"]),8)
            self.assertTrue(all(p["error_type"]=="SurfaceNotFoundError" for p in report["points"]))

    def test_parallel_equals_serial(self):
        densities=np.geomspace(1e17,1e19,7);eos=DegenerateNeutronGas()
        serial=sweep_central_densities(eos,densities)
        parallel=sweep_central_densities(eos,densities,workers=2)
        self.assertTrue(parallel.is_complete)
        for a,b in zip(serial.points,parallel.points):
            self.assertAlmostEqual(a.model.mass_kg/b.model.mass_kg,1.,places=14)
            self.assertAlmostEqual(a.model.radius_m/b.model.radius_m,1.,places=14)

    def test_optional_ui_export_without_streamlit_dependency(self):
        blob=result_zip([solve_star(DegenerateNeutronGas(),1e18)])
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            self.assertIn("star_1_summary.json",z.namelist())
            summary=json.loads(z.read("star_1_summary.json"))
            self.assertTrue(.62<summary["gravitational_mass_solar"]<.64)
            self.assertIn("manifest.json",z.namelist())


if __name__=="__main__":
    unittest.main()
