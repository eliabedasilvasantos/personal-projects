"""Quantified numerical validation, with analytic and published benchmarks.

Passing these checks verifies the stated mathematical models/regimes, not
the astrophysical realism of any uncalibrated EOS. No observations are bundled.
"""

from dataclasses import replace
import json
from math import pi, sqrt
from pathlib import Path
import numpy as np
from scipy.integrate import simpson

from .config import SolverConfig
from .constants import C, G
from .equations_of_state import DegenerateNeutronGas, IncompressibleEOS, RelativisticPolytrope, TabulatedEOS
from .stellar_properties import StellarModel, json_safe
from .tov_solver import solve_star, tov_rhs_si

OV_REFERENCE = {
    "title":"On Massive Neutron Cores", "authors":"J. R. Oppenheimer and G. M. Volkoff",
    "year":1939, "journal":"Physical Review 55, 374–381", "doi":"10.1103/PhysRev.55.374",
    "url":"https://doi.org/10.1103/PhysRev.55.374", "table":"Table I, finite t0 entries",
    "parameter_mapping":"x_c=sinh(t0/4), rho0_c=m_n^4 c^3 x_c^3/(3 pi^2 hbar^3)",
    "reported_precision":"historical rounded values, not high-precision reference data",
}


def dimensional_analysis() -> dict:
    """Exponent algebra in (kg,m,s), independent of numerical rescaling."""
    mass=np.array([1,0,0]); length=np.array([0,1,0]); time=np.array([0,0,1])
    gravity=-mass+3*length-2*time; speed=length-time
    energy_density=mass-length-2*time; pressure=energy_density.copy()
    pressure_mass=3*length+pressure-2*speed
    dm=2*length+energy_density-2*speed
    dp=gravity+energy_density+mass-2*speed-2*length
    passed=bool(np.array_equal(pressure_mass,mass) and np.array_equal(dm,mass-length) and
                np.array_equal(dp,pressure-length))
    return {"passed":passed,"basis":["kg","m","s"],"dm_dr":dm.tolist(),
            "dP_dr":dp.tolist(),"4pi_r3_P_c2":pressure_mass.tolist()}


def incompressible_benchmark() -> dict:
    rho,R=7e17,11000.
    eps=rho*C*C; mass=4*pi*rho*R**3/3; comp=G*mass/(R*C*C)
    a=sqrt(1-2*comp); pc=eps*(1-a)/(3*a-1)
    cfg=SolverConfig(rtol=1e-11,atol_mass=1e-14,atol_enthalpy=1e-14,
                     atol_baryonic_mass=1e-14,enforce_causality=False)
    model=solve_star(IncompressibleEOS(eps),config=cfg,central_pressure_pa=pc)
    r=model.solution.radius_m
    b=np.sqrt(1-2*comp*(r/R)**2)
    exact=eps*(b-a)/(3*a-b)
    errors={"mass_relative_error":abs(model.mass_kg/mass-1),
            "radius_relative_error":abs(model.radius_m/R-1),
            "pressure_profile_max_error_over_Pc":float(np.max(np.abs(model.solution.pressure_pa-exact))/pc)}
    return {"passed":all(v<1e-8 for v in errors.values()),"relative_acceptance":1e-8,
            "reference":"Schwarzschild interior constant-energy-density solution, 1916",
            "source":"https://arxiv.org/abs/physics/9912033",
            "status":"exact mathematical benchmark; acausal, not nuclear matter",
            "expected_mass_kg":mass,"expected_radius_m":R,**errors}


def newtonian_n1_benchmark() -> dict:
    # Analytic Lane-Emden n=1: theta=sin(xi)/xi; R=pi a; M=4 pi² a³ rho_c.
    K=.01; eos=RelativisticPolytrope(K,2)
    a=sqrt(K/(2*pi*G)); expected_R=pi*a
    rows=[]
    for sigma in (1e-4,1e-5,1e-6):
        rho=sigma*C*C/K
        cfg=SolverConfig(rtol=1e-11,atol_mass=1e-15,atol_enthalpy=1e-15,
                         atol_baryonic_mass=1e-15,max_step_m=60)
        model=solve_star(eos,rho,cfg)
        expected_M=4*pi*pi*a**3*rho
        xi=model.solution.radius_m/a
        theta=np.sinc(xi/np.pi)  # sin(xi)/xi with the correct center limit
        mask=model.solution.radius_m<.95*expected_R
        profile_error=float(np.max(np.abs(model.solution.rest_density_kg_m3[mask]/rho-theta[mask])))
        rows.append({"sigma_Pc_rho0c_c2":sigma,"mass_relative_error":abs(model.mass_kg/expected_M-1),
                     "radius_relative_error":abs(model.radius_m/expected_R-1),
                     "density_profile_error_over_rhoc":profile_error})
    passed=all(max(row[k] for k in ("mass_relative_error","radius_relative_error","density_profile_error_over_rhoc"))
               <20*row["sigma_Pc_rho0c_c2"] for row in rows)
    return {"passed":passed,"reference":"Analytic Newtonian Lane-Emden n=1 limit",
            "acceptance":"deviations <20 sigma; GR corrections should vanish as sigma→0",
            "expected_newtonian_radius_m":expected_R,"models":rows,
            "status":"asymptotic check, finite-sigma TOV is not exactly Newtonian"}


def fermi_published_benchmark() -> dict:
    eos=DegenerateNeutronGas(); rows=[]
    for t0,mass,radius in ((1,.30,21.1),(2,.60,13.3),(3,.71,9.5),(4,.64,6.8)):
        rho=eos.rho_scale*np.sinh(t0/4)**3
        model=solve_star(eos,float(rho),SolverConfig(rtol=1e-10,atol_mass=1e-13,
                          atol_enthalpy=1e-13,atol_baryonic_mass=1e-13,max_step_m=80))
        rows.append({"t0":t0,"central_rho0_kg_m3":rho,
            "published_mass_solar":mass,"computed_mass_solar":model.mass_solar,
            "published_radius_km":radius,"computed_radius_km":model.radius_km,
            "mass_relative_difference":abs(model.mass_solar/mass-1),
            "radius_relative_difference":abs(model.radius_km/radius-1)})
    return json_safe({"passed":all(max(r["mass_relative_difference"],r["radius_relative_difference"])<.04 for r in rows),
        "relative_acceptance":.04,"reference":OV_REFERENCE,"models":rows,
        "interpretation":"Percent-level reproduction of a rounded historical table using modern constants; not a 4% integrator error or observational validation."})


def profile_consistency(model: StellarModel) -> dict:
    """Independent quadrature plus fourth-order finite-difference TOV residuals.

    Differential residuals exclude r<0.02R and r>0.95R: division by vanishing
    gradients and nonsmooth surface powers would make relative errors ill-defined.
    Whole-star mass quadratures still include the center and the surface.
    """
    s=model.solution; r=s.radius_m; p=s.pressure_pa; m=s.gravitational_mass_kg; e=s.energy_density_j_m3
    integral=simpson(4*pi*r*r*e/C**2,x=r)
    # Samples after r0 are uniformly spaced; the origin-to-r0 interval is separate.
    dx=r[3]-r[2]
    dp=(p[1:-4]-8*p[2:-3]+8*p[4:-1]-p[5:])/(12*dx)
    dm=(m[1:-4]-8*m[2:-3]+8*m[4:-1]-m[5:])/(12*dx)
    idx=np.arange(3,len(r)-2)
    rhs=np.array([tov_rhs_si(float(r[i]),float(m[i]),float(p[i]),float(e[i])) for i in idx])
    mask=(r[idx]>.02*r[-1])&(r[idx]<.95*r[-1])
    dm_error=float(np.max(np.abs((dm[mask]-rhs[mask,0])/rhs[mask,0])))
    dp_error=float(np.max(np.abs((dp[mask]-rhs[mask,1])/rhs[mask,1])))
    mb_error=None
    if model.baryonic_mass_kg is not None:
        mb=simpson(4*pi*r*r*s.rest_density_kg_m3/np.sqrt(1-2*model.compactness_profile),x=r)
        mb_error=abs(mb/model.baryonic_mass_kg-1)
    mass_error=abs(integral/model.mass_kg-1)
    return json_safe({"passed":mass_error<2e-5 and dm_error<1e-4 and dp_error<1e-4 and
                                  (mb_error is None or mb_error<2e-5),
        "mass_integral_relative_error":mass_error,"baryonic_integral_relative_error":mb_error,
        "dm_dr_max_relative_residual":dm_error,"dP_dr_max_relative_residual":dp_error,
        "differential_domain":"0.02R<r<0.95R",
        "mass_integral_acceptance":2e-5,"differential_relative_acceptance":1e-4,
        "warning":"Finite-difference residuals depend on profile sampling, separate from integration tolerance."})


def _difference(model: StellarModel, reference: StellarModel) -> dict:
    return {"mass_solar":model.mass_solar,"radius_km":model.radius_km,
            "mass_relative_error":abs(model.mass_kg/reference.mass_kg-1),
            "radius_relative_error":abs(model.radius_m/reference.radius_m-1)}


def convergence_study() -> dict:
    eos=DegenerateNeutronGas();rho=1e18
    refcfg=SolverConfig(rtol=1e-12,atol_mass=1e-15,atol_enthalpy=1e-15,
                        atol_baryonic_mass=1e-15,max_step_m=20,profile_points=1401)
    reference=solve_star(eos,rho,refcfg);rows=[]
    for tol in (1e-5,1e-7,1e-9,1e-11):
        cfg=SolverConfig(rtol=tol,atol_mass=tol*1e-3,atol_enthalpy=tol*1e-3,
                         atol_baryonic_mass=tol*1e-3,max_step_m=500)
        model=solve_star(eos,rho,cfg)
        rows.append({"rtol":tol,"dimensionless_atol":tol*1e-3,
                     "function_evaluations":model.solution.diagnostics["function_evaluations"],**_difference(model,reference)})
    starts=[]
    for radius in (32.,8.,2.,.5,.125):
        model=solve_star(eos,rho,replace(refcfg,initial_radius_m=radius))
        starts.append({"initial_radius_m":radius,**_difference(model,reference)})
    scales=[]
    for L in (1e3,1e4,1e5):
        model=solve_star(eos,rho,replace(refcfg,length_scale_m=L))
        scales.append({"length_scale_m":L,**_difference(model,reference)})
    steps=[]
    for step in (1000.,300.,100.,30.):
        model=solve_star(eos,rho,replace(refcfg,max_step_m=step))
        steps.append({"max_step_m":step,**_difference(model,reference)})
    tight=max(rows[-1]["mass_relative_error"],rows[-1]["radius_relative_error"])
    loose=max(rows[0]["mass_relative_error"],rows[0]["radius_relative_error"])
    passed=tight<1e-7 and loose<2e-4 and all(max(r["mass_relative_error"],r["radius_relative_error"])<1e-6 for r in starts+scales+steps)
    return {"passed":passed,"central_rho0_kg_m3":rho,"reference":reference.summary(),
            "tolerance_study":{"models":rows,"tight_relative_acceptance":1e-7,"loose_relative_acceptance":2e-4},
            "center_study":starts,"unit_scale_study":scales,"max_step_study":steps,
            "other_studies_relative_acceptance":1e-6,
            "warning":"Differences to a refined numerical solution are not rigorous error bounds; strict monotonic decrease is not required near rounding/interpolation limits."}


def table_interpolation_study() -> dict:
    rho_c=1e18;K=.01; poly=RelativisticPolytrope(K,2)
    cfg=SolverConfig(rtol=1e-10,atol_mass=1e-13,atol_enthalpy=1e-13,atol_baryonic_mass=1e-13,max_step_m=60)
    reference=solve_star(poly,rho_c,cfg);source_rows=[];enthalpy_rows=[]
    for n in (41,81,161,321):
        rho=np.geomspace(1e10,1e19,n)
        table=TabulatedEOS(K*rho**2,rho*C*C+K*rho**2,rho,tail_gamma=2,
                            provenance="SYNTHETIC analytic-polytrope interpolation validation")
        model=solve_star(table,rho_c,cfg)
        source_rows.append({"source_rows":n,**_difference(model,reference),**table.interpolation_diagnostics})
    rho=np.geomspace(1e10,1e19,321)
    for n in (1025,4097,16385):
        table=TabulatedEOS(K*rho**2,rho*C*C+K*rho**2,rho,tail_gamma=2,
                 provenance="SYNTHETIC enthalpy-grid validation",enthalpy_grid_points=n)
        model=solve_star(table,rho_c,cfg)
        enthalpy_rows.append({"requested_enthalpy_nodes":n,**_difference(model,reference),**table.interpolation_diagnostics})
    passed=(max(source_rows[-1]["mass_relative_error"],source_rows[-1]["radius_relative_error"])<2e-5 and
            source_rows[-1]["mass_relative_error"]<source_rows[0]["mass_relative_error"] and
            enthalpy_rows[-1]["enthalpy_identity_max_relative_error"]<1e-4)
    return {"passed":passed,"reference":"Analytic Gamma=2 cold polytrope, synthetic table only",
            "source_grid_refinement":source_rows,"enthalpy_grid_refinement":enthalpy_rows,
            "fine_source_MR_relative_acceptance":2e-5,
            "warning":"This tests the importer, not any external nuclear EOS or a physical crust."}


def run_validation(output: str | Path | None = None) -> dict:
    """Run offline, deterministic physics/numerics checks and return a report."""
    checks={"dimensions":dimensional_analysis(),"incompressible":incompressible_benchmark(),
            "newtonian_n1":newtonian_n1_benchmark(),"published_fermi":fermi_published_benchmark(),
            "convergence":convergence_study(),"table_interpolation":table_interpolation_study(),
            "profile_consistency":profile_consistency(solve_star(DegenerateNeutronGas(),1e18))}
    report={"passed":all(c["passed"] for c in checks.values()),"checks":checks,
        "scope":"Validated mathematical benchmark models and numerical consistency in the tested regimes.",
        "not_validated":["Realistic nuclear matter, composition and observational inference",
            "Radial-mode eigenvalues and nonlinear time evolution", "Rotation, tides, magnetic fields, finite temperature",
            "Imported scientific EOS tables (none redistributed), density-jump interfaces"],
        "observational_data_included":False}
    if output is not None:
        path=Path(output);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(json_safe(report),indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    return report
