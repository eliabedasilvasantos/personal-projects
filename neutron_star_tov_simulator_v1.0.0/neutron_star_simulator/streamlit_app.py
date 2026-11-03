"""Optional UI; the scientific core is independent of Streamlit.

Run from the package's parent directory:
  streamlit run neutron_star_simulator/streamlit_app.py
"""

import io
import json
from pathlib import Path
import sys
import tempfile
import zipfile

# Streamlit can execute this file directly without installing the package.
if __package__ in {None,""}:
    sys.path.insert(0,str(Path(__file__).resolve().parent.parent))

import numpy as np

from neutron_star_simulator.constants import C, M_NEUTRON
from neutron_star_simulator.config import SolverConfig
from neutron_star_simulator.equations_of_state import (
    DegenerateNeutronGas,RelativisticPolytrope,PiecewisePolytrope,ConstantSoundSpeed,TabulatedEOS,
)
from neutron_star_simulator.exceptions import SimulatorError
from neutron_star_simulator.stability import sweep_central_densities,analyze_stability,refine_maximum,export_sequence
from neutron_star_simulator.tov_solver import solve_star
from neutron_star_simulator.main import write_manifest
from neutron_star_simulator.visualization import plot_profiles,plot_sequences


def result_zip(models) -> bytes:
    """Downloadable reproducible summaries and profiles for compared stars."""
    buffer=io.BytesIO()
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        for i,model in enumerate(models): model.export(root,f"star_{i+1}")
        write_manifest(root,["streamlit","star comparison"],extra={"models":len(models)})
        with zipfile.ZipFile(buffer,"w",compression=zipfile.ZIP_DEFLATED) as archive:
            for file in sorted(root.iterdir()): archive.write(file,file.name)
    return buffer.getvalue()


def render_app() -> None:
    try:
        import streamlit as st
    except ImportError as exc:
        raise SystemExit("Install optional UI dependencies: pip install -r neutron_star_simulator/requirements-ui.txt") from exc
    st.set_page_config(page_title="TOV • Estrutura relativística",layout="wide")
    st.title("Estrutura relativística de estrelas compactas")
    st.caption("Equilíbrios frios e esféricos • massa gravitacional, raio, perfis e sequências TOV")
    st.info("Resultados calculados de modelos. A estabilidade radial não é demonstrada: os modos de pulsação não são resolvidos.")
    st.sidebar.header("Equação de estado")
    choice=st.sidebar.selectbox("EOS",["Fermi de nêutrons livres","Politrópica fria","Politrópica por trechos","Velocidade do som constante","Tabela CSV"])
    eos=None
    try:
        if choice=="Fermi de nêutrons livres":
            eos=DegenerateNeutronGas()
            st.sidebar.caption("T=0, g=2; sem interações nucleares, equilíbrio beta ou crosta.")
        elif choice=="Politrópica fria":
            st.sidebar.caption("Valores iniciais são parâmetros sintéticos de experimento, sem ajuste nuclear.")
            gamma=st.sidebar.number_input("Gamma",min_value=1.01,max_value=5.,value=2.,step=.05)
            lr=st.sidebar.number_input("log10 rho_ref [kg/m³]",value=17.,step=.1)
            lp=st.sidebar.number_input("log10 P_ref [Pa]",value=32.,step=.1)
            eos=RelativisticPolytrope.from_reference(10**lr,10**lp,gamma)
            st.sidebar.caption(f"K = {eos.K:.6g} em unidades SI dependentes de Gamma")
        elif choice=="Politrópica por trechos":
            st.sidebar.caption("Exemplo sintético contínuo; não é uma EOS nuclear publicada.")
            k=st.sidebar.number_input("K0 em SI",value=float(.01*(1e17)**.4),format="%.6e")
            gammas=[float(v) for v in st.sidebar.text_input("Gammas","1.6,2.2,1.8").split(",")]
            breaks=[float(v) for v in st.sidebar.text_input("Quebras rho0 [kg/m³]","1e17,8e17").split(",")]
            eos=PiecewisePolytrope(k,gammas,breaks)
        elif choice=="Velocidade do som constante":
            st.sidebar.caption("Matéria autoligada parametrizada; composição não determinada.")
            sound=st.sidebar.number_input("c_s²/c²",min_value=.01,max_value=1.,value=1/3,step=.01)
            le=st.sidebar.number_input("log10 epsilon_s [J/m³]",value=float(np.log10(4e34)))
            lr=st.sidebar.number_input("log10 rho0_s [kg/m³]",value=float(np.log10(4e34/C**2)))
            eos=ConstantSoundSpeed(sound,10**le,10**lr)
        else:
            uploaded=st.sidebar.file_uploader("CSV canônico frio 1D",type="csv")
            units=st.sidebar.selectbox("Unidades",["si","cgs","nuclear"])
            tail=st.sidebar.number_input("Gamma da extensão até P=0",min_value=1.01,value=1.5)
            provenance=st.sidebar.text_area("Fonte, composição, temperatura, convenção de energia")
            baryon=st.sidebar.number_input("Massa bariônica [kg] (para unidades nucleares)",value=M_NEUTRON,format="%.12e")
            if uploaded and provenance.strip():
                with tempfile.TemporaryDirectory() as tmp:
                    path=Path(tmp)/"uploaded_eos.csv";path.write_bytes(uploaded.getvalue())
                    eos=TabulatedEOS.from_csv(path,units=units,tail_gamma=tail,provenance=provenance,baryon_mass_kg=baryon)
                st.sidebar.json(eos.interpolation_diagnostics)
            else:
                st.warning("Selecione a tabela e informe sua proveniência para executar.")
        with st.sidebar.expander("Controle numérico"):
            lr_tol=st.number_input("log10 rtol",min_value=-12.,max_value=-3.,value=-8.)
            la_tol=st.number_input("log10 atol adimensional",min_value=-16.,max_value=-5.,value=-11.)
            r0=st.number_input("Raio inicial [m]",min_value=.001,value=.5)
            rmax=st.number_input("Raio máximo [km]",min_value=1.,value=100.)
            step=st.number_input("Passo máximo [m]",min_value=1.,value=150.)
            method=st.selectbox("Método",["DOP853","RK45","Radau"])
        cfg=SolverConfig(method=method,rtol=10**lr_tol,atol_mass=10**la_tol,atol_enthalpy=10**la_tol,
                         atol_baryonic_mass=10**la_tol,initial_radius_m=r0,max_radius_m=rmax*1000,max_step_m=step)
    except (SimulatorError,ValueError,OverflowError) as exc:
        st.error(f"EOS/configuração inválida: {exc}");return
    if eos is None: return
    st.session_state.setdefault("models",[])
    single,sequence_tab=st.tabs(["Estrela e comparação","Sequência massa–raio"])
    with single:
        logc=st.number_input("log10 densidade central de REPOUSO [kg/m³]",min_value=8.,max_value=22.,value=18.,step=.05)
        if st.button("Integrar estrela",type="primary"):
            try:
                with st.spinner("Integrando TOV e localizando a superfície..."):
                    model=solve_star(eos,10**logc,cfg)
                st.session_state.models.append(model)
                st.session_state.models=st.session_state.models[-6:]
            except (SimulatorError,ValueError,OverflowError) as exc: st.error(f"Integração rejeitada: {exc}")
        models=st.session_state.models
        if models:
            model=models[-1]
            cols=st.columns(4)
            for col,label,value in zip(cols,["M [M☉]","R [km]","GM/(Rc²)","Redshift z"],
                                       [model.mass_solar,model.radius_km,model.compactness,model.surface_redshift]):
                col.metric(label,f"{value:.6g}")
            st.dataframe([{"EOS":m.solution.eos_metadata["name"],"rho0,c [kg/m³]":m.solution.central_rest_density_kg_m3,
                           "M [M_sun]":m.mass_solar,"R [km]":m.radius_km,"C":m.compactness} for m in models])
            with tempfile.TemporaryDirectory() as tmp:
                figures=plot_profiles(models,tmp)
                st.image(Path([f for f in figures if f.endswith("profiles_comparison.png")][0]).read_bytes())
            st.download_button("Baixar perfis e metadados ZIP",result_zip(models),"tov_models.zip",mime="application/zip")
            with st.expander("Metadados da última estrela"): st.json(model.summary())
            if st.button("Limpar comparação"):
                st.session_state.models=[];st.rerun()
    with sequence_tab:
        c1,c2,c3=st.columns(3)
        lo=c1.number_input("log10 rho0 mínima",value=16.)
        hi=c2.number_input("log10 rho0 máxima",value=19.3)
        count=c3.number_input("Número de estrelas",min_value=8,max_value=250,value=50,step=1)
        assumed=st.checkbox("Assumo que o início da sequência pertence a um ramo radialmente estável",value=False)
        if st.button("Calcular sequência"):
            try:
                if hi<=lo: raise ValueError("A densidade máxima deve exceder a mínima")
                with st.spinner("Varredura adaptativa em andamento..."):
                    sequence=sweep_central_densities(eos,np.geomspace(10**lo,10**hi,int(count)),cfg)
                    peak=refine_maximum(eos,sequence)
                st.session_state["sequence_result"]=(sequence,peak,assumed)
            except (SimulatorError,ValueError,OverflowError) as exc: st.error(str(exc))
        if "sequence_result" in st.session_state:
            sequence,peak,assumed=st.session_state["sequence_result"]
            failures=[{"rho0":p.requested_rest_density_kg_m3,"erro":p.error} for p in sequence.points if p.model is None]
            if failures: st.warning(f"{len(failures)} pontos falharam e permanecem registrados.");st.dataframe(failures)
            with tempfile.TemporaryDirectory() as tmp:
                files=plot_sequences({sequence.eos_metadata["name"]:sequence},tmp,assume_low_density_stable=assumed)
                st.image(Path([f for f in files if f.endswith("mass_radius.png")][0]).read_bytes())
                export_sequence(sequence,tmp,assume_low_density_stable=assumed,maximum=peak)
                write_manifest(Path(tmp),["streamlit","density sequence"])
                buffer=io.BytesIO()
                with zipfile.ZipFile(buffer,"w",zipfile.ZIP_DEFLATED) as z:
                    for path in sorted(Path(tmp).iterdir()): z.write(path,path.name)
                st.download_button("Baixar sequência e gráficos",buffer.getvalue(),"tov_sequence.zip",mime="application/zip")
            st.json(peak)
            with st.expander("Critério de estabilidade e limitações"):
                st.json(analyze_stability(sequence,assume_low_density_stable=assumed))


if __name__=="__main__":
    render_app()
