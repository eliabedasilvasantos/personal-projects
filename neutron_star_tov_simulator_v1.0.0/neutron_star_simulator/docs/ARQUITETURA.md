# Arquitetura e execução por etapas

O diretório neutron_star_simulator é simultaneamente pacote Python e raiz
da distribuição. Execute os comandos da pasta que contém esse diretório.

| Módulo | Responsabilidade |
|---|---|
| constants.py | Constantes SI, conversões, escala geométrica |
| config.py | Parâmetros numéricos validados e serializáveis |
| equations_of_state.py | EOSState, interface EOS, EOS analíticas e tabeladas |
| tov_solver.py | TOV SI, TOV adimensional, série central e evento superficial |
| stellar_properties.py | Perfil imutável por convenção, propriedades e exportações |
| stability.py | Varredura, falhas por ponto, derivadas e refino de máximos |
| visualization.py | Figuras científicas PNG/PDF sem dependência da interface |
| validation.py | Benchmarks analíticos, literatura, convergência e resíduos |
| main.py | CLI reprodutível, resultados, hashes e manifesto de ambiente |
| streamlit_app.py | Interface opcional que reutiliza o núcleo |
| exceptions.py | Exceções físicas e numéricas explícitas |
| tests/ | Testes executáveis com unittest, sem dados remotos |
| data/ | Esquema e tabela sintética para testar importação |
| results/ | Saídas reais de execuções de validação e demonstração |

Fluxo: entrada SI → EOS e configuração → integrador → StellarModel →
validação/propriedades → sequência e estabilidade → figuras/CSV/JSON.
O núcleo não depende de Streamlit, Numba ou JAX. Varreduras independentes
podem usar processos; não há estado global mutável nem números aleatórios.

Etapa 1: formulação, unidades e configuração; testar conversões.
Etapa 2: EOS; testar limites, inversões, primeira lei, continuidade e causalidade.
Etapa 3: integrador; testar centro, superfície e solução exata.
Etapa 4: propriedades; testar quadraturas e identidades integrais.
Etapa 5: sequência; testar ordenação, derivadas, máximos e falhas.
Etapa 6: gráficos; gerar e inspecionar figuras com dados calculados.
Etapa 7: validação; executar literatura, convergência e regressões físicas.
Etapa 8: interface e paralelismo; testar CLI, importação opcional e igualdade
serial/paralela. O registro final está em docs/VALIDACAO.md e nos JSONs.
