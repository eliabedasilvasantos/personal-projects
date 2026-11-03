"""Derived observables and portable results, separate from integration."""

from dataclasses import dataclass
import json
from math import isfinite, sqrt
from pathlib import Path
from typing import Any
import numpy as np

from . import __version__
from .constants import C, G, M_SUN_KG
from .tov_solver import TOVSolution


def json_safe(value: Any) -> Any:
    """Strict JSON: unavailable quantities become null, never NaN/Infinity."""
    if isinstance(value, dict):
        return {str(k):json_safe(v) for k,v in value.items()}
    if isinstance(value, (list,tuple,np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float,np.floating)):
        return float(value) if isfinite(value) else None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


@dataclass(frozen=True)
class StellarModel:
    solution: TOVSolution

    @property
    def radius_m(self) -> float:
        return float(self.solution.radius_m[-1])

    @property
    def radius_km(self) -> float:
        return self.radius_m/1000

    @property
    def mass_kg(self) -> float:
        """Total mass = gravitational mass; includes energy and gravitational binding."""
        return float(self.solution.gravitational_mass_kg[-1])

    @property
    def mass_solar(self) -> float:
        return self.mass_kg/M_SUN_KG

    @property
    def baryonic_mass_kg(self) -> float | None:
        mass = float(self.solution.baryonic_mass_kg[-1])
        return mass if isfinite(mass) else None

    @property
    def binding_energy_j(self) -> float | None:
        return None if self.baryonic_mass_kg is None else (self.baryonic_mass_kg-self.mass_kg)*C*C

    @property
    def compactness(self) -> float:
        """C=GM/(R c²); 2C is explicitly a separate quantity."""
        return G*self.mass_kg/(self.radius_m*C*C)

    @property
    def surface_redshift(self) -> float:
        return 1/sqrt(1-2*self.compactness)-1

    @property
    def compactness_profile(self) -> np.ndarray:
        s = self.solution
        return np.divide(G*s.gravitational_mass_kg,C*C*s.radius_m,
                         out=np.zeros_like(s.radius_m),where=s.radius_m>0)

    @property
    def lapse_profile(self) -> np.ndarray:
        """exp(Phi) from dPhi=-dh, matched to Schwarzschild vacuum at R."""
        return sqrt(1-2*self.compactness)*np.exp(-self.solution.enthalpy)

    def summary(self) -> dict:
        s = self.solution
        return json_safe({
            "software_version":__version__, "eos":s.eos_metadata,
            "radius_m":self.radius_m,"radius_km":self.radius_km,
            "gravitational_mass_kg":self.mass_kg,"gravitational_mass_solar":self.mass_solar,
            "total_mass_definition":"total mass is the gravitational mass m(R)",
            "baryonic_mass_kg":self.baryonic_mass_kg,
            "binding_energy_j":self.binding_energy_j,
            "compactness_GM_Rc2":self.compactness,
            "schwarzschild_ratio_2GM_Rc2":2*self.compactness,
            "surface_redshift":self.surface_redshift,
            "central_rest_density_kg_m3":s.central_rest_density_kg_m3,
            "central_energy_density_j_m3":s.central_energy_density_j_m3,
            "central_mass_energy_density_kg_m3":s.central_energy_density_j_m3/C**2,
            "central_pressure_pa":s.central_pressure_pa,"central_enthalpy":s.central_enthalpy,
            "configuration":s.config.to_dict(),"diagnostics":s.diagnostics,
            "interpretation":"Computed equilibrium, not an observation; radial modes not solved.",
        })

    def export(self, directory: str | Path, stem: str = "star") -> dict[str,str]:
        directory = Path(directory)
        directory.mkdir(parents=True,exist_ok=True)
        if Path(stem).name != stem or stem in {"", ".", ".."}:
            raise ValueError("stem must be a simple filename prefix")
        s = self.solution
        data = np.c_[s.radius_m,s.gravitational_mass_kg,s.pressure_pa,s.energy_density_j_m3,
                     s.rest_density_kg_m3,s.baryonic_mass_kg,s.enthalpy,self.compactness_profile,
                     s.sound_speed_squared_over_c2,self.lapse_profile]
        csv_path, json_path = directory/f"{stem}_profile.csv", directory/f"{stem}_summary.json"
        header = ("radius_m,gravitational_mass_kg,pressure_pa,epsilon_j_m3,rho0_kg_m3,"
                  "baryonic_mass_kg,enthalpy,compactness_Gm_rc2,cs2_over_c2,lapse")
        np.savetxt(csv_path,data,delimiter=",",header=header,comments="",fmt="%.16e")
        json_path.write_text(json.dumps(self.summary(),indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
        return {"profile_csv":str(csv_path),"summary_json":str(json_path)}
