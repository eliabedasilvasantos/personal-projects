"""Publication-style Matplotlib figures, with SI-labelled physical data.

All figures are computed equilibria, never observational constraints. PNG and
vector PDF are saved together. This module does not change global rcParams.
"""

from pathlib import Path
from typing import Mapping, Sequence
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from .constants import C, M_SUN_KG, MEV_FM3_TO_PA
from .equations_of_state import EquationOfState
from .exceptions import SimulatorError
from .stability import StellarSequence, analyze_stability
from .stellar_properties import StellarModel

COLORS = ("#1f5a85", "#bb5a16", "#39735b", "#824b8b", "#5f6268", "#a62539")
STYLE = {"font.family":"serif","font.serif":["STIXGeneral","DejaVu Serif"],
         "mathtext.fontset":"stix","font.size":11,"axes.labelsize":12,
         "axes.titlesize":13,"legend.fontsize":9,"axes.spines.top":False,
         "axes.spines.right":False,"axes.linewidth":.8,"xtick.direction":"in",
         "ytick.direction":"in","savefig.dpi":240,"figure.dpi":110,
         "lines.linewidth":1.8,"pdf.fonttype":42,"ps.fonttype":42}


def _save(fig, prefix: Path) -> list[str]:
    prefix.parent.mkdir(parents=True,exist_ok=True)
    paths=[]
    for suffix in (".png", ".pdf"):
        path=prefix.with_suffix(suffix)
        fig.savefig(path,bbox_inches="tight",metadata={"Creator":"neutron_star_simulator"} if suffix==".pdf" else None)
        paths.append(str(path))
    plt.close(fig)
    return paths


def _finish(ax, *, xlabel: str, ylabel: str, title: str, logy: bool=False):
    ax.set(xlabel=xlabel,ylabel=ylabel,title=title)
    if logy: ax.set_yscale("log")
    ax.grid(True,which="major",alpha=.18,linewidth=.6)
    ax.legend(frameon=False)


def plot_profiles(models: Sequence[StellarModel], output: str | Path) -> list[str]:
    """Compare central densities for one EOS or explicitly labelled EOSs."""
    if not models: raise ValueError("At least one stellar model is required")
    output=Path(output); paths=[]
    definitions=(
        ("density",r"$\rho_0,\;\varepsilon/c^2$ [kg m$^{-3}$]","Densidade de repouso e de massa-energia",True),
        ("pressure",r"$P$ [Pa]","Pressão relativística",True),
        ("enclosed_mass",r"$m(r)$ [$M_\odot$]","Massa gravitacional acumulada",False),
        ("compactness_profile",r"$Gm(r)/(rc^2)$","Compactação no interior",False),
    )
    with plt.rc_context(STYLE):
        for key,ylabel,title,logy in definitions:
            fig,ax=plt.subplots(figsize=(6.8,4.6),layout="constrained")
            for i,model in enumerate(models):
                s=model.solution; r=s.radius_m/1000; color=COLORS[i%len(COLORS)]
                label=rf"$\rho_{{0,c}}={s.central_rest_density_kg_m3:.2e}$ kg m$^{{-3}}$"
                if key=="density":
                    keep=s.rest_density_kg_m3>0
                    ax.plot(r[keep],s.rest_density_kg_m3[keep],color=color,label=label)
                    keep=s.energy_density_j_m3>0
                    ax.plot(r[keep],s.energy_density_j_m3[keep]/C**2,color=color,ls="--",alpha=.7)
                elif key=="pressure":
                    keep=s.pressure_pa>0
                    ax.plot(r[keep],s.pressure_pa[keep],color=color,label=label)
                elif key=="enclosed_mass":
                    ax.plot(r,s.gravitational_mass_kg/M_SUN_KG,color=color,label=label)
                else:
                    ax.plot(r,model.compactness_profile,color=color,label=label)
            _finish(ax,xlabel=r"Raio areal $r$ [km]",ylabel=ylabel,title=title,logy=logy)
            if key=="density":
                ax.text(.02,.03,r"linha cheia: $\rho_0$; tracejada: $\varepsilon/c^2$",
                        transform=ax.transAxes,fontsize=9)
            fig.suptitle("Equilíbrios TOV calculados • sem dados observacionais",fontsize=10,color="#555555")
            paths += _save(fig,output/key)
        fig,axes=plt.subplots(2,2,figsize=(10,7.5),layout="constrained")
        for i,model in enumerate(models):
            s=model.solution;r=s.radius_m/1000;color=COLORS[i%len(COLORS)]
            label=rf"$\rho_{{0,c}}={s.central_rest_density_kg_m3:.2e}$"
            keep=s.rest_density_kg_m3>0
            axes[0,0].semilogy(r[keep],s.rest_density_kg_m3[keep],color=color,label=label)
            keep=s.pressure_pa>0
            axes[0,1].semilogy(r[keep],s.pressure_pa[keep],color=color,label=label)
            axes[1,0].plot(r,s.gravitational_mass_kg/M_SUN_KG,color=color,label=label)
            axes[1,1].plot(r,model.compactness_profile,color=color,label=label)
        for ax,ylabel,title in zip(axes.flat,[r"$\rho_0$ [kg m$^{-3}$]",r"$P$ [Pa]",r"$m$ [$M_\odot$]",r"$Gm/(rc^2)$"],
                                   ["Densidade de repouso","Pressão","Massa acumulada","Compactação"]):
            _finish(ax,xlabel="r [km]",ylabel=ylabel,title=title)
        fig.suptitle("Comparação de modelos e densidades centrais",fontsize=15)
        paths += _save(fig,output/"profiles_comparison")
    return paths


def plot_sequences(sequences: Mapping[str,StellarSequence], output: str | Path, *,
                   assume_low_density_stable: bool=False) -> list[str]:
    output=Path(output); paths=[]
    with plt.rc_context(STYLE):
        for mode in ("mass_radius","mass_central_density","compactness_sequence"):
            fig,ax=plt.subplots(figsize=(7.5,5.2),layout="constrained")
            for i,(label,sequence) in enumerate(sequences.items()):
                color=COLORS[i%len(COLORS)]
                analysis=analyze_stability(sequence,assume_low_density_stable=assume_low_density_stable)
                labelled=False
                for j in range(len(sequence.points)-1):
                    first,second=sequence.points[j].model,sequence.points[j+1].model
                    if first is None or second is None: continue  # visible gaps, no bridging failures
                    pair=[first,second]
                    if mode=="mass_radius":
                        xx=[m.radius_km for m in pair]; yy=[m.mass_solar for m in pair]
                    elif mode=="mass_central_density":
                        xx=[m.solution.central_energy_density_j_m3/C**2 for m in pair]; yy=[m.mass_solar for m in pair]
                    else:
                        xx=[m.mass_solar for m in pair]; yy=[m.compactness for m in pair]
                    state=analysis["labels"][j]
                    style="--" if "negative_slope" in state else ":" if "unresolved" in state else "-"
                    ax.plot(xx,yy,color=color,ls=style,label=label if not labelled else None)
                    labelled=True
                peak=analysis["largest_sample_index"]
                if peak is not None and analysis["largest_sample_is_bracketed_peak"]:
                    model=sequence.points[peak].model
                    if mode=="mass_radius": xx,yy=model.radius_km,model.mass_solar
                    elif mode=="mass_central_density": xx,yy=model.solution.central_energy_density_j_m3/C**2,model.mass_solar
                    else: xx,yy=model.mass_solar,model.compactness
                    ax.plot(xx,yy,"o",color=color,ms=4)
            if mode=="mass_radius":
                xlabel,ylabel,title="Raio R [km]",r"Massa gravitacional M [$M_\odot$]","Relação massa–raio"
            elif mode=="mass_central_density":
                xlabel,ylabel,title=r"$\varepsilon_c/c^2$ [kg m$^{-3}$]",r"$M$ [$M_\odot$]","Sequências e pontos de retorno"
                ax.set_xscale("log")
            else:
                xlabel,ylabel,title=r"$M$ [$M_\odot$]",r"$GM/(Rc^2)$","Compactação superficial"
            _finish(ax,xlabel=xlabel,ylabel=ylabel,title=title)
            fig.suptitle("EOS de referência / parametrizações não calibradas",fontsize=11,color="#555555")
            fig.text(.11,-.035,"Tracejado: derivada negativa; pontilhado: extremo/inconclusivo.\n"
                     "Os estilos não substituem uma análise dos modos radiais.",fontsize=8)
            paths += _save(fig,output/mode)
    return paths


def plot_equations_of_state(equations: Mapping[str,EquationOfState], output: str | Path) -> list[str]:
    with plt.rc_context(STYLE):
        fig,axes=plt.subplots(1,2,figsize=(10.3,4.4),layout="constrained")
        for i,(label,eos) in enumerate(equations.items()):
            data=[]
            for rho in np.geomspace(1e16,2e19,250):
                try:
                    s=eos.state_from_rest_density(float(rho))
                    data.append((s.energy_density_j_m3/C**2,s.pressure_pa/MEV_FM3_TO_PA,s.sound_speed_squared_over_c2))
                except (SimulatorError,ValueError): continue
            if not data: continue
            a=np.array(data);color=COLORS[i%len(COLORS)]
            axes[0].loglog(a[:,0],a[:,1],label=label,color=color)
            axes[1].semilogx(a[:,0],a[:,2],label=label,color=color)
        axes[1].axhline(1,color="black",ls="--",lw=1,label="limite causal")
        _finish(axes[0],xlabel=r"$\varepsilon/c^2$ [kg m$^{-3}$]",ylabel=r"$P$ [MeV fm$^{-3}$]",title="Pressão versus energia")
        _finish(axes[1],xlabel=r"$\varepsilon/c^2$ [kg m$^{-3}$]",ylabel=r"$c_s^2/c^2$",title="Compressibilidade e causalidade")
        fig.suptitle("EOS analíticas e parametrizadas • parâmetros de experimento",fontsize=13)
        return _save(fig,Path(output)/"equations_of_state")


def plot_convergence(report: dict, output: str | Path) -> list[str]:
    rows=report["tolerance_study"]["models"]
    with plt.rc_context(STYLE):
        fig,ax=plt.subplots(figsize=(6.8,4.5),layout="constrained")
        x=np.array([r["rtol"] for r in rows]);
        for key,label,color in (("mass_relative_error","M",COLORS[0]),("radius_relative_error","R",COLORS[1])):
            y=np.array([max(r[key],1e-16) for r in rows])
            ax.loglog(x,y,"o-",color=color,label=label)
        _finish(ax,xlabel="Tolerância relativa do integrador",ylabel="Diferença relativa à integração de referência",title="Sensibilidade numérica — gás livre")
        ax.invert_xaxis()
        ax.text(.03,.03,"Valores abaixo de $10^{-16}$ limitados apenas para exibição.\n"
                "Refino sem monotonicidade estrita pode atingir arredondamento.",transform=ax.transAxes,fontsize=8)
        return _save(fig,Path(output)/"convergence")
