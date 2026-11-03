# Simulador relativístico de estrelas de nêutrons — TOV

Projeto Python executável para equilíbrios estáticos de estrelas compactas
frias, com EOS intercambiáveis, superfície por evento, perfis internos,
sequências massa–raio, diagnóstico de pontos de retorno e validação quantitativa.

**O que está validado:** implementação dos modelos matemáticos nos regimes
testados, incluindo solução exata, limite newtoniano e referência histórica.
**O que não está validado:** microfísica nuclear realista, inferência
observacional ou estabilidade dinâmica por modos de pulsação. Os parâmetros
dos exemplos politrópicos/CSS são sintéticos, não ajustes a dados.

## 1. Começar

Requer Python >=3.10. A execução entregue usou Python 3.12.14, NumPy 2.3.5,
SciPy 1.17.0 e Matplotlib 3.10.8 em precisão float64. Não há dados remotos
necessários para os testes ou exemplos. As faixas de versões permitidas não
significam que todas as combinações foram testadas.

Depois de extrair o ZIP, abra o terminal **na pasta que contém o diretório
neutron_star_simulator**. Todos os comandos desta seção partem dessa pasta.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r neutron_star_simulator/requirements.txt
python -m unittest discover -s neutron_star_simulator/tests -v
python -m neutron_star_simulator.main validate --out results/validation
```

Para reproduzir as versões exatas, use requirements-tested.txt em lugar de
requirements.txt. No Windows, a ativação é `.venv\Scripts\activate`.
Não execute simplesmente `python neutron_star_simulator/main.py`: a CLI usa
imports de pacote. Use a forma `python -m ...` acima.

Uma estrela de referência, com **densidade central de repouso** de 1e18 kg/m³:

```bash
python -m neutron_star_simulator.main star --eos fermi --rho-c 1e18 --out results/fermi_star
```

Resultado esperado nas versões testadas: M≈0.628616 M_sun, R≈12.670080 km,
C≈0.0733. Pequenas diferenças de arredondamento são normais. Os valores são
de um gás livre de nêutrons, não de uma EOS nuclear realista.

Varredura e refino do máximo:

```bash
python -m neutron_star_simulator.main sweep --eos fermi --rho-min 1e16 --rho-max 2e19 --points 60 --out results/fermi_sequence
```

Comparações e todas as figuras, inteiramente offline após instalar bibliotecas:

```bash
python -m neutron_star_simulator.main demo --points 60 --out results/demo
```

Paralelismo por processos, com resultados na mesma ordem da malha:

```bash
python -m neutron_star_simulator.main sweep --eos fermi --rho-min 1e16 --rho-max 2e19 --points 100 --workers 2 --out results/parallel
```

Para chamadas paralelas em scripts próprios, use a guarda
`if __name__ == "__main__":`. Em notebooks/ambientes interativos onde spawn
não consegue importar o módulo principal, use workers=1 ou execute a CLI.

Instalação opcional como pacote, da mesma pasta:

```bash
python -m pip install ./neutron_star_simulator
neutron-star-simulator --help
```

## 2. Fundamento matemático

Leia primeiro [docs/FUNDAMENTACAO.md](docs/FUNDAMENTACAO.md), depois
[docs/EOS.md](docs/EOS.md) e [docs/ARQUITETURA.md](docs/ARQUITETURA.md).
São parte do projeto, não apenas comentários no código.

O problema usa as equações completas:

\[
\frac{dm}{dr}=4\pi r^2\varepsilon/c^2,\qquad
\frac{dP}{dr}=-\frac{G(\varepsilon+P)(m+4\pi r^3P/c^2)}
{c^2r^2(1-2Gm/(rc^2))}.
\]

rho0 é densidade de repouso (kg/m³), epsilon é energia TOTAL por volume
(J/m³) e P é pressão (Pa). Não se confundem rho0 e epsilon/c². Uma barótropa
fria fornece a relação entre as três, com as hipóteses de composição declaradas.

Para resolver a superfície, usa-se a transformação exata
h=integral[dP/(epsilon+P)] e dh/dr=-(G/c²)(m+4pi r³P/c²)/[r²(1-2Gm/(rc²))].
As variáveis integradas são adimensionais (mu,h,mu_b); solve_ivp/DOP853
usa controle adaptativo e evento terminal h=0. A origem é tratada por série
regular, nunca por divisão em r=0. TOV em SI também está exposta e testada
contra a forma adimensional.

## 3. EOS disponíveis

| Modelo | Conteúdo e limites |
|---|---|
| DegenerateNeutronGas | Gás de Fermi relativístico exato em T=0, sem interações; referência Oppenheimer–Volkoff |
| RelativisticPolytrope | P=K rho0^Gamma, epsilon=rho0 c²+P/(Gamma-1), primeira lei fria |
| PiecewisePolytrope | Pressão, energia e potencial químico contínuos; sem transições com salto |
| ConstantSoundSpeed | P=s(epsilon-epsilon_s), matéria autoligada parametrizada; normalização bariônica explícita |
| TabulatedEOS | CSV SI/CGS/nuclear, PCHIP monótono, domínio finito, proveniência e extensão superficial declarada |
| IncompressibleEOS | Somente teste acausal contra solução exata de Schwarzschild |

Politrópica de experimento (K=0.01 SI, Gamma=2, **sem calibração nuclear**):

```bash
python -m neutron_star_simulator.main star --eos polytrope --k 0.01 --gamma 2 --rho-c 1e18 --out results/poly
```

Parâmetros podem ser gravados em JSON, usando os campos de
data/polytrope_experiment.json:

```bash
python -m neutron_star_simulator.main star --eos-json neutron_star_simulator/data/polytrope_experiment.json --rho-c 1e18 --out results/poly_json
```

O importador funciona sem tabelas científicas proprietárias/anexadas. Uma
tabela **sintética** politrópica está incluída somente para exercitar a rota:

```bash
python -m neutron_star_simulator.main star --eos table --table-path neutron_star_simulator/data/synthetic_polytrope_si.csv --table-units si --tail-gamma 2 --provenance "SYNTHETIC Gamma=2 K=0.01 SI; not nuclear data" --rho-c 1e18 --out results/table_test
```

Para fontes científicas reais, seleção, download, conversão e auditoria,
leia [data/README.md](data/README.md). Nenhuma tabela nuclear externa foi
validada nesta entrega; não há EOS SLy/APR/DD2 implementada sob esses nomes.

## 4. Propriedades e arquivos gerados

Uma estrela exporta star_profile.csv e star_summary.json, com:

- Raio areal e massa gravitacional total M=m(R).
- Perfis de m, P, epsilon, rho0, entalpia e velocidade do som.
- Massa bariônica em volume próprio e energia de ligação, quando definidas.
- C=GM/(Rc²), 2GM/(Rc²), redshift superficial e lapso exp(Phi).
- Condições centrais, parâmetros EOS, tolerâncias, evento e diagnósticos.

Massa total e massa gravitacional são a mesma grandeza neste problema;
massa bariônica é apresentada separadamente. A unidade solar usa o parâmetro
nominal IAU dividido por G.

Uma sequência exporta sequence.csv e sequence.json, inclusive pontos que
falharam e mensagens de erro. Pontos ausentes interrompem curvas e derivadas.
O máximo é refinado com novas integrações; se estiver na borda da malha ou
ao lado de uma falha, não se afirma ter encontrado a massa máxima.

Cada comando gera manifest.json com versões, argumentos e hashes SHA-256
do código. Tabelas importadas também registram seu hash e proveniência.
Arquivos numéricos CSV conservam 16 dígitos de escrita, que não representam
16 dígitos de precisão física. JSON não usa NaN/Infinity; valores indisponíveis
aparecem como null. Em CSV, quantidades não definidas ou falhas usam NaN.

Os comandos sobrescrevem os arquivos de mesmo nome no diretório escolhido.
Use `--out` diferente para preservar experimentos anteriores.

## 5. Estabilidade e gráficos

São gerados PNG e PDF vetorial para densidade, pressão, massa acumulada,
compactação radial e superficial, massa–raio, massa–densidade central,
comparação entre EOS, velocidades do som e convergência.

As curvas logarítmicas omitem P=0 e rho0=0 somente na representação gráfica;
a raiz superficial permanece nos dados exportados. Cores/estilos diferenciam
modelos e diagnósticos, não dados observacionais. Os pontos nos gráficos M–R
marcam máximos **amostrados**; os valores refinados estão nos JSONs.

Por padrão uma derivada positiva é classificada como “não classificada”.
Somente use `--assume-low-density-stable` quando houver fundamento para
essa hipótese no início da sequência. Mesmo assim a classificação é de
**candidato** no ramo inicial, não uma prova dinâmica. Derivadas negativas,
regiões de extremo e ramos crescentes posteriores têm rótulos próprios.
Não se resolvem autovalores radiais nem se pressupõe que uma terceira família
seja estável apenas por sua inclinação positiva.

## 6. Controle numérico e falhas

As opções --rtol e --atol atuam sobre variáveis adimensionais. Também estão
expostos --max-step-m, --initial-radius-m, --max-radius-km, --length-scale-m,
--profile-points e --method. Aumentar profile-points melhora a amostragem
e as diferenças finitas, não a trajetória integrada. Reduzir rtol não corrige
EOS mal resolvida, convenções erradas ou microfísica omitida.

Causalidade c_s<=c é exigida em todas as avaliações por padrão. O override
--allow-acausal existe para experimentos explícitos e fica nos metadados.
Estrelas sem evento superficial ou com geometria imprópria não são retornadas
como soluções. A interface mantém o último resultado válido separado de um
erro novo; os metadados correspondem ao resultado que foi efetivamente calculado.

Saídas da CLI: 0=êxito dos controles aplicáveis, 1=erro ou falha de validação,
2=sequência incompleta ou diagnóstico de perfil fora dos limites de aceitação.
Resíduos de diferenças finitas perto de interfaces EOS podem exigir amostragem
específica; consulte profile_checks.json, sem ocultar o problema.

## 7. Interface opcional

```bash
python -m pip install -r neutron_star_simulator/requirements-ui.txt
streamlit run neutron_star_simulator/streamlit_app.py
```

Permite selecionar EOS, alterar densidade e parâmetros, ajustar tolerâncias,
integrar estrelas, guardar até seis modelos para comparação, fazer varreduras,
ver falhas e baixar resultados. Os gráficos mostram o resultado após cada
integração adaptativa, sem animação que sugira evolução temporal física.

**Estado da verificação:** módulo e exportação ZIP foram testados sem exigir
Streamlit. O Streamlit não está instalado no ambiente da entrega; a navegação
visual da interface não foi executada. Instalação e ensaio interativo devem
ser feitos no ambiente de uso. O núcleo, a CLI, o paralelismo e as figuras
foram efetivamente executados.

## 8. Uso pela API

```python
import numpy as np
from neutron_star_simulator.config import SolverConfig
from neutron_star_simulator.equations_of_state import DegenerateNeutronGas
from neutron_star_simulator.tov_solver import solve_star
from neutron_star_simulator.stability import sweep_central_densities, analyze_stability

eos = DegenerateNeutronGas()
star = solve_star(eos, central_rest_density_kg_m3=1e18,
                  config=SolverConfig(rtol=1e-9, atol_mass=1e-12,
                                      atol_enthalpy=1e-12, atol_baryonic_mass=1e-12))
print(star.mass_solar, star.radius_km, star.surface_redshift)
star.export("results/api")
sequence = sweep_central_densities(eos, np.geomspace(1e16, 2e19, 60))
diagnostic = analyze_stability(sequence, assume_low_density_stable=False)
```

Para EOS incompressível, a condição central deve ser pressão, pois a energia
constante não determina P_c. Na API, forneça central_pressure_pa e configuração
com enforce_causality=False. É somente um benchmark matemático acausal.

## 9. Documentação, testes e expansões

Leia [docs/VALIDACAO.md](docs/VALIDACAO.md) para resultados efetivos e
limitações, [docs/PERFORMANCE.md](docs/PERFORMANCE.md) para Numba/JAX e
paralelismo, [docs/REFERENCIAS.md](docs/REFERENCIAS.md) para fontes primárias
verificadas e [docs/ETAPAS.md](docs/ETAPAS.md) para a sequência de implementação.

Extensões físicas futuras exigem equações e validações adicionais: crosta
nuclear acoplada, transições com salto, autovalores de pulsações radiais,
deformabilidade tidal/Love, rotação lenta, composição e temperatura. Não há
funções vazias ou resultados fictícios apresentados como essas extensões.
