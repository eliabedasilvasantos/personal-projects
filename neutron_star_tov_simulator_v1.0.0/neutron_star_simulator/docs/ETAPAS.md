# Etapas de desenvolvimento e estado atual

## Etapa 1 — Fundamentação matemática e física — concluída

As equações TOV completas, a métrica, as hipóteses de fluido perfeito, a
separação entre rho0 e epsilon, a adimensionalização e a expansão regular do
centro estão em `FUNDAMENTACAO.md`. Os testes de conversão SI e invariância de
escala estão em `tests/test_units.py` e `tests/test_solver.py`.

## Etapa 2 — Equações de estado — concluída

Foram implementados gás livre de Fermi, politrópica fria, politrópica por
trechos, CSS e importador CSV estrito. A primeira lei fria, continuidade,
limites de Fermi, som e conversões são testados em `tests/test_eos.py`.
Parâmetros de exemplos são sintéticos e explicitamente não calibrados.

## Etapa 3 — Integrador TOV — concluída

`solve_ivp` com DOP853 por padrão, tolerância por componente, evento de
superfície h=0, evento de horizonte, diagnóstico de passos e método alternativo
foram implementados em `tov_solver.py`. A solução analítica incompressível é
reproduzida com erro inferior a 1e-8.

## Etapa 4 — Propriedades estelares — concluída

`stellar_properties.py` calcula M, R, massa bariônica, energia de ligação,
compactação, redshift, lapso e exportação CSV/JSON. Quadraturas independentes
e serialização sem NaN/Infinity são testadas.

## Etapa 5 — Varredura e estabilidade — concluída com escopo limitado

`stability.py` constrói sequências, mantém falhas, calcula derivadas, localiza
extremos amostrados e refina máximos bracketed. A classificação é cautelosa;
modos radiais não são resolvidos. O próximo avanço científico seria implementar
o problema de autovalor das pulsações radiais e tratar interfaces com salto.

## Etapa 6 — Visualização — concluída

PNG e PDF vetorial para perfis, massa–raio, densidade central, compactação,
EOS e convergência são gerados sem depender de observações. O teste de figuras
confirma cabeçalhos PNG/PDF e tamanho não vazio.

## Etapa 7 — Testes e validação — concluída para os regimes declarados

49 testes automatizados passam no ambiente da entrega. O relatório inclui
dimensionalidade, solução exata, limite newtoniano, comparação histórica,
convergência, tabelas sintéticas e resíduos de perfil. Consulte
`VALIDACAO.md` para fronteiras da inferência.

## Etapa 8 — Interface e otimização — núcleo concluído; UI opcional

CLI, exportação ZIP, paralelismo por processos e o módulo Streamlit estão
implementados. O núcleo e o exportador da UI foram testados; Streamlit não
está instalado no ambiente atual, portanto a navegação visual da página não
foi executada. Numba e JAX permanecem extensões experimentais documentadas em
`PERFORMANCE.md`, sem dependência ou ganho de velocidade alegado.
