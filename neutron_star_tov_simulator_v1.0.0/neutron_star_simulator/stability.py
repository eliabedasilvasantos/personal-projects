"""Central-density sequences and cautious turning-point diagnostics.

This module never calls a positive mass slope a proof of dynamic stability.
Failures remain in place; derivatives are never taken across missing models.
"""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any
import multiprocessing
import numpy as np
from scipy.optimize import minimize_scalar

from .config import SolverConfig
from .constants import M_SUN_KG
from .equations_of_state import EquationOfState
from .exceptions import SimulatorError
from .stellar_properties import StellarModel, json_safe
from .tov_solver import solve_star


@dataclass(frozen=True)
class SweepPoint:
    requested_rest_density_kg_m3: float
    model: StellarModel | None
    error_type: str | None = None
    error: str | None = None


@dataclass(frozen=True)
class StellarSequence:
    points: tuple[SweepPoint, ...]
    eos_metadata: dict[str, Any]

    @property
    def successful_models(self) -> list[StellarModel]:
        return [p.model for p in self.points if p.model is not None]

    @property
    def is_complete(self) -> bool:
        return all(p.model is not None for p in self.points)


def _one_point(args) -> SweepPoint:
    eos, rho, config = args
    try:
        return SweepPoint(rho, solve_star(eos, rho, config))
    except (SimulatorError, ValueError, OverflowError, FloatingPointError) as exc:
        return SweepPoint(rho, None, type(exc).__name__, str(exc))


def sweep_central_densities(eos: EquationOfState, densities_kg_m3,
                           config: SolverConfig | None = None, *, workers: int = 1) -> StellarSequence:
    densities = np.asarray(densities_kg_m3, dtype=float)
    if (densities.ndim != 1 or len(densities)<1 or not np.all(np.isfinite(densities)) or
        np.any(densities<=0) or np.any(np.diff(densities)<=0)):
        raise ValueError("Central densities must be a finite, strictly increasing positive 1D grid")
    if not isinstance(workers,int) or workers<1:
        raise ValueError("workers must be a positive integer")
    args = [(eos,float(rho),config or SolverConfig()) for rho in densities]
    if workers==1:
        points = tuple(_one_point(a) for a in args)
    else:
        # spawn avoids inheriting BLAS/GUI state; callers need the __main__ guard.
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context("spawn")) as pool:
            points = tuple(pool.map(_one_point,args))
    return StellarSequence(points,eos.metadata())


def analyze_stability(sequence: StellarSequence, *, assume_low_density_stable: bool = False,
                      slope_tolerance: float = 1e-3) -> dict:
    """Compute dM/d(epsilon_c) on contiguous segments of >=3 equilibria.

    Candidate initial stability is assigned only under an EXPLICIT supplied
    assumption about the low-density starting branch. Finite grids can miss
    extrema; no mode counting or first-order-transition criterion is included.
    """
    if not 0 < slope_tolerance < 1:
        raise ValueError("slope_tolerance must be in (0,1)")
    n = len(sequence.points)
    labels = np.full(n,"failed",dtype=object)
    derivatives = np.full(n,np.nan)
    log_slopes = np.full(n,np.nan)
    valid = np.array([p.model is not None for p in sequence.points])
    indices = np.flatnonzero(valid)
    groups = np.split(indices,np.where(np.diff(indices)>1)[0]+1) if len(indices) else []
    extrema: list[dict] = []
    for group in groups:
        if not len(group): continue
        labels[group] = "insufficient_contiguous_models"
        if len(group)<3: continue
        models = [sequence.points[int(i)].model for i in group]
        eps = np.array([m.solution.central_energy_density_j_m3 for m in models])
        mass = np.array([m.mass_kg for m in models])
        if np.any(np.diff(eps)<=0):
            raise ValueError("Central energy density must rise along a fixed-EOS sequence")
        dlog = np.gradient(mass,np.log(eps),edge_order=2)
        derivatives[group] = dlog/eps
        log_slopes[group] = dlog/max(mass)
        local_extrema = []
        for j in range(1,len(group)-1):
            kind = ("maximum" if mass[j]>mass[j-1] and mass[j]>mass[j+1] else
                    "minimum" if mass[j]<mass[j-1] and mass[j]<mass[j+1] else None)
            if kind:
                item = {"index":int(group[j]),"type":kind,
                        "central_rest_density_kg_m3":sequence.points[int(group[j])].requested_rest_density_kg_m3,
                        "mass_solar":models[j].mass_solar,
                        "bracket_indices":[int(group[j-1]),int(group[j+1])]}
                extrema.append(item); local_extrema.append(j)
        first = min(local_extrema) if local_extrema else len(group)
        for j, idx in enumerate(group):
            slope = log_slopes[idx]
            if abs(slope)<slope_tolerance or j in local_extrema:
                labels[idx] = "turning_region_unresolved"
            elif slope<0:
                labels[idx] = "negative_slope_candidate_unstable"
            elif j>first:
                labels[idx] = "rising_after_turning_unresolved"
            elif group[0]==0 and assume_low_density_stable:
                labels[idx] = "initial_branch_candidate_stable"
            else:
                labels[idx] = "positive_slope_unclassified"
    successful = [i for i,p in enumerate(sequence.points) if p.model]
    largest = max(successful,key=lambda i:sequence.points[i].model.mass_kg) if successful else None
    bracketed = (largest is not None and 0<largest<n-1 and
                 all(sequence.points[j].model is not None for j in (largest-1,largest+1)) and
                 sequence.points[largest].model.mass_kg > sequence.points[largest-1].model.mass_kg and
                 sequence.points[largest].model.mass_kg > sequence.points[largest+1].model.mass_kg)
    return json_safe({"labels":labels.tolist(),"dM_depsilon_c_kg_m3_per_j":derivatives,
        "dM_dlogepsilon_over_maxM":log_slopes,"sampled_extrema":extrema,
        "largest_sample_index":largest,"largest_sample_is_bracketed_peak":bool(bracketed),
        "assume_low_density_stable":assume_low_density_stable,
        "slope_tolerance":slope_tolerance,"complete_sequence":sequence.is_complete,
        "limitations":["No radial pulsation eigenvalues were computed.",
            "An initial stable branch is an assumption, not an output of the derivative test.",
            "A rising branch after extrema needs mode counting; positive slope alone is insufficient.",
            "A finite interval/grid does not establish a global maximum or resolve all extrema.",
            "No derivative is computed across a failed model."]})


def refine_maximum(eos: EquationOfState, sequence: StellarSequence,
                   config: SolverConfig | None = None, *, index: int | None = None,
                   log_density_tolerance: float = 1e-5) -> dict:
    """Refine a bracketed sampled local maximum by NEW TOV integrations.

    If the largest model lies on a boundary or beside a failure, return a
    censored status, never extrapolate a maximum. Multiple maxima are refined
    individually by supplying index; the highest sampled mass is the default.
    """
    if eos.metadata() != sequence.eos_metadata:
        raise ValueError("Peak refinement must use the same EOS as the sequence")
    if not 0 < log_density_tolerance < 1:
        raise ValueError("log_density_tolerance must be in (0,1)")
    valid = [i for i,p in enumerate(sequence.points) if p.model]
    if not valid:
        return {"status":"no_successful_models","model":None}
    idx = index if index is not None else max(valid,key=lambda i:sequence.points[i].model.mass_kg)
    if idx not in valid:
        raise ValueError("index must identify a successful model")
    if idx==0 or idx==len(sequence.points)-1:
        return {"status":"boundary_limited_no_maximum_claim","model":sequence.points[idx].model.summary()}
    a,b,c = sequence.points[idx-1:idx+2]
    if a.model is None or c.model is None:
        return {"status":"failure_limited_no_maximum_claim","model":b.model.summary()}
    if not b.model.mass_kg>max(a.model.mass_kg,c.model.mass_kg):
        return {"status":"not_a_bracketed_maximum","model":b.model.summary()}
    cfg = config or b.model.solution.config
    bounds = np.log([a.requested_rest_density_kg_m3,c.requested_rest_density_kg_m3])
    def objective(logrho):
        return -solve_star(eos,float(np.exp(logrho)),cfg).mass_solar
    try:
        opt = minimize_scalar(objective,bounds=tuple(bounds),method="bounded",
                              options={"xatol":log_density_tolerance,"maxiter":100})
        if not opt.success:
            return {"status":"optimizer_failed","message":str(opt.message),"model":None}
        model = solve_star(eos,float(np.exp(opt.x)),cfg)
    except (SimulatorError, ValueError, OverflowError) as exc:
        return {"status":"refinement_integration_failed","message":str(exc),"model":None}
    if not model.mass_kg>=b.model.mass_kg*(1-10*cfg.rtol):
        return {"status":"refinement_inconsistent_with_sample","model":model.summary()}
    return {"status":"refined_local_maximum","sample_index":idx,
            "density_bracket_kg_m3":[a.requested_rest_density_kg_m3,c.requested_rest_density_kg_m3],
            "log_density_optimizer_tolerance":log_density_tolerance,
            "optimizer_evaluations":int(opt.nfev),"model":model.summary(),
            "interpretation":"Local peak inside the supplied interval; no dynamic stability proof."}


def export_sequence(sequence: StellarSequence, directory: str | Path, stem: str="sequence", *,
                    assume_low_density_stable: bool=False, maximum: dict | None=None) -> dict:
    directory = Path(directory); directory.mkdir(parents=True,exist_ok=True)
    if Path(stem).name!=stem or stem in {"", ".", ".."}:
        raise ValueError("stem must be a simple filename prefix")
    analysis = analyze_stability(sequence,assume_low_density_stable=assume_low_density_stable)
    # CSV retains every request; failed rows use NaN and explicit status in JSON.
    rows=[]
    for i,p in enumerate(sequence.points):
        m=p.model
        rows.append([p.requested_rest_density_kg_m3,
            m.solution.central_energy_density_j_m3 if m else np.nan,
            m.mass_solar if m else np.nan, m.radius_km if m else np.nan,
            m.compactness if m else np.nan,m.surface_redshift if m else np.nan,
            analysis["dM_depsilon_c_kg_m3_per_j"][i] if m and analysis["dM_depsilon_c_kg_m3_per_j"][i] is not None else np.nan])
    csv_path=directory/f"{stem}.csv"; json_path=directory/f"{stem}.json"
    np.savetxt(csv_path,rows,delimiter=",",fmt="%.16e",comments="",
        header="central_rho0_kg_m3,central_epsilon_j_m3,mass_solar,radius_km,compactness_GM_Rc2,surface_redshift,dM_depsilon_c_kg_m3_per_j")
    report={"eos":sequence.eos_metadata,"stability":analysis,"maximum":maximum,
            "points":[{"requested_rest_density_kg_m3":p.requested_rest_density_kg_m3,
                       "model":p.model.summary() if p.model else None,
                       "error_type":p.error_type,"error":p.error} for p in sequence.points]}
    json_path.write_text(json.dumps(json_safe(report),indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    return {"sequence_csv":str(csv_path),"sequence_json":str(json_path)}
