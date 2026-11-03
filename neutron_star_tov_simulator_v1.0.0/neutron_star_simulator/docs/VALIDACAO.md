# Registro de validação

Os comandos abaixo foram executados no ambiente da entrega em 2 de novembro
de 2026. O relatório completo, com números e tolerâncias, está em
`results/validation/validation.json` no pacote final. O registro não confunde
erro numérico com incerteza da física.

## Resultado executado

| Verificação | Resultado | O que ela cobre |
|---|---:|---|
| Dimensões | PASS | Potências SI das duas equações TOV |
| Esfera incompressível | PASS | Solução interior exata de Schwarzschild; erro de massa, raio e perfil < 1e-8 |
| Limite newtoniano n=1 | PASS | EOS Gamma=2 em três valores decrescentes de sigma |
| Tabela Oppenheimer–Volkoff | PASS | Quatro casos finitos da Tabela I, diferença arredondada < 4% |
| Convergência | PASS | rtol, centro inicial, escala geométrica e passo máximo |
| Interpolação tabulada | PASS | Densidade da tabela, primeira lei e malha de entalpia |
| Perfil independente | PASS | Quadraturas de massa/bárions e resíduos finitos de dM/dr, dP/dr |
| Testes automatizados | PASS | 49 testes `unittest`, incluindo EOS, TOV, sequências, CLI e figuras |

O modelo livre de nêutrons reproduziu aproximadamente M_max=0,71018 M_sun e
R=9,16 km no ponto de retorno refinado. Isso é consistente com a referência
histórica quando se usa a EOS de gás livre, mas não descreve matéria nuclear
realista. A fonte original informa valores arredondados; a diferença de poucos
por cento é uma comparação de referência publicada, não o erro do integrador.

## O que foi verificado numericamente

As grandezas derivadas foram recalculadas por quadratura de Simpson a partir
dos perfis, independentemente do estado de massa integrado. A massa bariônica
foi recalculada com o fator de volume próprio
`(1-2Gm/(rc²))^-1/2`. Os resíduos diferenciais usam diferenças de quarta ordem
apenas em `0,02R < r < 0,95R`, evitando divisão por gradientes que desaparecem
no centro/superfície.

A convergência compara contra uma integração refinada. Ela mostra estabilidade
da resposta às tolerâncias, ao passo máximo, ao raio inicial da expansão e à
escala geométrica. Uma mudança menor que a tolerância do integrador não é um
limite de erro rigoroso. EOS tabuladas têm ainda erro de fonte e interpolação,
que ficam registrados separadamente.

## Limites da validação

O registro não prova que uma EOS sintética representa a matéria nuclear. Não há
EOS CompOSE, SLy, APR, DD2 ou outra tabela nuclear externa redistribuída ou
calibrada nesta versão. Não foram calculados autovalores de pulsação radial,
estabilidade dinâmica, rotação, deformabilidade tidal, campos magnéticos,
temperatura finita, composição beta ou interfaces de primeira ordem.

O critério dM/d epsilon_c só rotula candidatos sob as hipóteses documentadas
em `FUNDAMENTACAO.md`; uma inclinação positiva depois de um extremo permanece
inconclusiva. A grade de densidades é finita e um pico na borda não recebe
status de massa máxima. Falhas de integração não são removidas das sequências.

## Reproduzir

```bash
python -m unittest discover -s neutron_star_simulator/tests -v
python -m neutron_star_simulator.main validate --out results/validation
python -m neutron_star_simulator.main demo --points 60 --out results/demo
```

Cada execução produz um `manifest.json` com versões das dependências,
argumentos e hashes SHA-256 dos módulos Python. Para uma publicação, fixe o
ambiente, preserve o manifesto, forneça a EOS e cite a fonte da tabela usada.
