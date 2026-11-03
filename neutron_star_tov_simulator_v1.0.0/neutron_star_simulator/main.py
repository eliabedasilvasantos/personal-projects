"""Command line experiments. Run: python -m neutron_star_simulator.main --help."""

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
from hashlib import sha256
import json
import platform
from pathlib import Path
import sys
import numpy as np
import scipy
import matplotlib

from . import __version__
from .config import SolverConfig
from .constants import C
from .equations_of_state import (
    EquationOfState, DegenerateNeutronGas, RelativisticPolytrope, PiecewisePolytrope,
    ConstantSoundSpeed, TabulatedEOS, eos_from_spec,
)
from .exceptions import SimulatorError
from .stability import sweep_central_densities, refine_maximum, export_sequence
from .stellar_properties import json_safe
from .tov_solver import solve_star
from .validation import run_validation, profile_consistency


def _write_json(path: Path, value: dict | list) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(json_safe(value),indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")


def write_manifest(directory: Path, arguments: list[str], *, extra: dict | None=None) -> Path:
    root=Path(__file__).parent
    hashes={str(p.relative_to(root)):sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*.py")) if "__pycache__" not in p.parts}
    manifest={"software_version":__version__,"created_at_utc":datetime.now(timezone.utc).isoformat(),
        "arguments":arguments,"python":sys.version,"platform":platform.platform(),
        "dependencies":{"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":matplotlib.__version__},
        "source_sha256":hashes,"extra":extra or {},
        "data_status":"Computed equilibria; no observational dataset is used."}
    path=directory/"manifest.json";_write_json(path,manifest);return path


def _eos_arguments(parser):
    parser.add_argument("--eos",choices=["fermi","polytrope","piecewise","css","table"],default="fermi")
    parser.add_argument("--eos-json",type=Path,help="Explicit EOS specification; overrides EOS CLI options")
    parser.add_argument("--k",type=float,help="Polytropic K in SI; for piecewise, K of the first segment")
    parser.add_argument("--gamma",type=float,help="Exponent of a cold polytrope")
    parser.add_argument("--gammas",type=str,help="Comma-separated piecewise exponents")
    parser.add_argument("--breaks",type=str,help="Comma-separated rest-density breaks, kg/m³")
    parser.add_argument("--sound-speed-squared",type=float,help="CSS: c_s²/c²")
    parser.add_argument("--epsilon-surface",type=float,help="CSS: surface energy density, J/m³")
    parser.add_argument("--rho0-surface",type=float,help="CSS: surface baryonic rest density, kg/m³")
    parser.add_argument("--table-path",type=Path)
    parser.add_argument("--table-units",choices=["si","cgs","nuclear"])
    parser.add_argument("--tail-gamma",type=float,help="Explicit low-density table extension")
    parser.add_argument("--provenance",help="Scientific source, conditions, energy convention and composition")
    parser.add_argument("--baryon-mass-kg",type=float,help="Mandatory for table nuclear units")


def _numeric_arguments(parser):
    parser.add_argument("--method",default="DOP853",choices=["DOP853","RK45","Radau","BDF","LSODA"])
    parser.add_argument("--rtol",type=float,default=1e-8)
    parser.add_argument("--atol",type=float,default=1e-11,help="Dimensionless atol for mu,h,mu_b")
    parser.add_argument("--initial-radius-m",type=float,default=.5)
    parser.add_argument("--max-radius-km",type=float,default=100.)
    parser.add_argument("--max-step-m",type=float,default=150.)
    parser.add_argument("--length-scale-m",type=float,default=10000.)
    parser.add_argument("--profile-points",type=int,default=701)
    parser.add_argument("--allow-acausal",action="store_true",help="EXPLICIT experimental override, recorded in output")


def _get_eos(args) -> EquationOfState:
    if args.eos_json:
        return eos_from_spec(json.loads(args.eos_json.read_text(encoding="utf-8")))
    def need(*keys):
        if any(getattr(args,key) is None for key in keys):
            raise ValueError("Explicit EOS parameters required: "+", ".join("--"+k.replace("_","-") for k in keys))
    if args.eos=="fermi": return DegenerateNeutronGas()
    if args.eos=="polytrope":
        need("k","gamma");return RelativisticPolytrope(args.k,args.gamma)
    if args.eos=="piecewise":
        need("k","gammas","breaks")
        return PiecewisePolytrope(args.k,[float(v) for v in args.gammas.split(",")],
                                  [float(v) for v in args.breaks.split(",")])
    if args.eos=="css":
        need("sound_speed_squared","epsilon_surface","rho0_surface")
        return ConstantSoundSpeed(args.sound_speed_squared,args.epsilon_surface,args.rho0_surface)
    need("table_path","table_units","tail_gamma","provenance")
    return TabulatedEOS.from_csv(args.table_path,units=args.table_units,tail_gamma=args.tail_gamma,
                         provenance=args.provenance,baryon_mass_kg=args.baryon_mass_kg)


def _get_config(args) -> SolverConfig:
    return SolverConfig(method=args.method,rtol=args.rtol,atol_mass=args.atol,
        atol_enthalpy=args.atol,atol_baryonic_mass=args.atol,
        initial_radius_m=args.initial_radius_m,max_radius_m=args.max_radius_km*1000,
        max_step_m=args.max_step_m,length_scale_m=args.length_scale_m,
        profile_points=args.profile_points,enforce_causality=not args.allow_acausal)


def demo(output: Path, *, points: int=60, workers: int=1) -> dict:
    """Offline benchmark experiments, not calibrated nuclear predictions."""
    from .visualization import plot_profiles, plot_sequences, plot_equations_of_state
    if points<8: raise ValueError("demo needs at least eight density samples")
    specs={
        "fermi":{"type":"fermi"},
        "poly_K005":{"type":"polytrope","K_SI":.005,"gamma":2.},
        "poly_K010":{"type":"polytrope","K_SI":.01,"gamma":2.},
        "poly_K020":{"type":"polytrope","K_SI":.02,"gamma":2.},
        "piecewise":{"type":"piecewise","K0_SI":.01*(1e17)**.4,
                     "gammas":[1.6,2.2,1.8],"breaks_kg_m3":[1e17,8e17]},
        "css":{"type":"css","s":1/3,"epsilon_surface_j_m3":4e34,"rho0_surface_kg_m3":4e34/C**2},
    }
    labels={"fermi":"Fermi livre","poly_K005":r"Politrópica $\Gamma=2$, $K=0.005$",
        "poly_K010":r"Politrópica $\Gamma=2$, $K=0.010$","poly_K020":r"Politrópica $\Gamma=2$, $K=0.020$",
        "piecewise":"Politrópica por trechos (sintética)","css":r"CSS $s=1/3$ (sintética)"}
    _write_json(output/"experimental_eos.json",{"status":"UNCALIBRATED experiment parameters; not published nuclear EOS fits",
        "K_units":"Pa/(kg m^-3)^Gamma; identical units among Gamma=2 examples",
        "specifications":specs})
    cfg=SolverConfig(profile_points=701)
    sequences={};equations={};peaks={}
    for key,spec in specs.items():
        eos=eos_from_spec(spec)
        lower=(eos.rho0_surface_kg_m3*1.01 if isinstance(eos,ConstantSoundSpeed) else 1e16)
        sequence=sweep_central_densities(eos,np.geomspace(lower,2e19,points),cfg,workers=workers)
        peak=refine_maximum(eos,sequence)
        # Low-density stability is NOT assumed for this comparison plot.
        export_sequence(sequence,output/key,maximum=peak)
        sequences[labels[key]]=sequence; equations[labels[key]]=eos;peaks[key]=peak
        failures=sum(p.model is None for p in sequence.points)
        print(f"{key}: {len(sequence.points)-failures}/{len(sequence.points)} equilibria; {peak['status']}")
    fermi=DegenerateNeutronGas()
    models=[solve_star(fermi,rho,cfg) for rho in (3e17,1e18,4e18)]
    for i,model in enumerate(models): model.export(output/"profiles",f"fermi_{i+1}")
    figures=plot_profiles(models,output/"figures")
    figures+=plot_sequences(sequences,output/"figures")
    figures+=plot_equations_of_state(equations,output/"figures")
    _write_json(output/"peak_comparison.json",peaks)
    return {"models_per_sequence":points,"sequences":len(specs),"figures":figures,
        "scope":"Reference EOS and synthetic uncalibrated parameters; no observational constraints."}


def build_parser() -> argparse.ArgumentParser:
    parser=argparse.ArgumentParser(description="Cold spherical relativistic stars: TOV in SI/geometrized units")
    sub=parser.add_subparsers(dest="command",required=True)
    star=sub.add_parser("star",help="Integrate one equilibrium and export profile, observables, plots")
    _eos_arguments(star);_numeric_arguments(star)
    central=star.add_mutually_exclusive_group(required=True)
    central.add_argument("--rho-c",type=float,help="Central REST density, kg/m³")
    central.add_argument("--central-pressure",type=float,help="Central pressure, Pa")
    star.add_argument("--out",type=Path,default=Path("results/star"))
    star.add_argument("--no-plots",action="store_true")
    sweep=sub.add_parser("sweep",help="Construct a central-density sequence and refine its sampled peak")
    _eos_arguments(sweep);_numeric_arguments(sweep)
    sweep.add_argument("--rho-min",type=float,required=True);sweep.add_argument("--rho-max",type=float,required=True)
    sweep.add_argument("--points",type=int,default=60);sweep.add_argument("--workers",type=int,default=1)
    sweep.add_argument("--assume-low-density-stable",action="store_true")
    sweep.add_argument("--out",type=Path,default=Path("results/sweep"));sweep.add_argument("--no-plots",action="store_true")
    validation=sub.add_parser("validate",help="Run offline analytic, published and numerical validation")
    validation.add_argument("--out",type=Path,default=Path("results/validation"))
    d=sub.add_parser("demo",help="Reproduce synthetic EOS comparisons and all profile figures")
    d.add_argument("--points",type=int,default=60);d.add_argument("--workers",type=int,default=1)
    d.add_argument("--out",type=Path,default=Path("results/demo"))
    return parser


def main(argv: list[str] | None=None) -> int:
    supplied=list(argv) if argv is not None else sys.argv[1:]
    args=build_parser().parse_args(supplied)
    output=args.out;output.mkdir(parents=True,exist_ok=True)
    try:
        if args.command=="validate":
            report=run_validation(output/"validation.json")
            from .visualization import plot_convergence
            plot_convergence(report["checks"]["convergence"],output/"figures")
            write_manifest(output,supplied,extra={"passed":report["passed"]})
            for name,check in report["checks"].items(): print(f"{name}: {'PASS' if check['passed'] else 'FAIL'}")
            return 0 if report["passed"] else 1
        if args.command=="demo":
            result=demo(output,points=args.points,workers=args.workers)
            write_manifest(output,supplied,extra=result);return 0
        eos=_get_eos(args);cfg=_get_config(args)
        if args.command=="star":
            star=solve_star(eos,args.rho_c,cfg,central_pressure_pa=args.central_pressure)
            star.export(output)
            checks=profile_consistency(star)
            _write_json(output/"profile_checks.json",checks)
            if not args.no_plots:
                from .visualization import plot_profiles
                plot_profiles([star],output/"figures")
            write_manifest(output,supplied,extra={"eos":eos.metadata(),"config":cfg.to_dict()})
            print(f"M = {star.mass_solar:.9f} M_sun; R = {star.radius_km:.9f} km; C = {star.compactness:.8f}; z = {star.surface_redshift:.8f}")
            if not checks["passed"]:
                print("Profile diagnostics exceeded thresholds; inspect profile_checks.json and refine sampling near EOS interfaces. Exit code 2.",file=sys.stderr)
                return 2
            return 0
        if args.points<3 or not 0<args.rho_min<args.rho_max:
            raise ValueError("sweep needs >=3 samples and 0<rho-min<rho-max")
        sequence=sweep_central_densities(eos,np.geomspace(args.rho_min,args.rho_max,args.points),cfg,workers=args.workers)
        peak=refine_maximum(eos,sequence)
        export_sequence(sequence,output,assume_low_density_stable=args.assume_low_density_stable,maximum=peak)
        if not args.no_plots:
            from .visualization import plot_sequences
            plot_sequences({eos.name:sequence},output/"figures",assume_low_density_stable=args.assume_low_density_stable)
        write_manifest(output,supplied,extra={"eos":eos.metadata(),"config":cfg.to_dict(),
                                            "complete_sequence":sequence.is_complete,"maximum_status":peak["status"]})
        print(f"{len(sequence.successful_models)}/{len(sequence.points)} equilibria; peak: {peak['status']}")
        if not sequence.is_complete:
            print("Sequence has explicit failures; inspect sequence.json. Exit code 2.",file=sys.stderr)
            return 2
        return 0
    except (SimulatorError,ValueError,KeyError,OverflowError,OSError) as exc:
        print(f"ERROR [{type(exc).__name__}]: {exc}",file=sys.stderr)
        return 1


if __name__=="__main__":
    raise SystemExit(main())
