# Equações de estado e fronteiras físicas

## Politrópica relativística fria

P=K rho0^Gamma; epsilon=rho0 c²+P/(Gamma-1), Gamma>1. A expressão para epsilon
é obtida da primeira lei fria d(epsilon/rho0)=-P d(1/rho0), com energia de
repouso c² no limite rho0→0. Portanto P=K(epsilon/c²)^Gamma seria uma EOS
diferente. Não se trocam as duas definições.

h=ln[1+Gamma K rho0^(Gamma-1)/((Gamma-1)c²)] e
c_s²/c²=dP/depsilon=Gamma P/(epsilon+P). Gamma>2 pode tornar a EOS acausal
em alta densidade; o integrador rejeita os estados infratores por padrão.
K tem unidades que dependem de Gamma. Comparações de rigidez utilizam
P_ref e rho_ref comuns, declarados no arquivo de experimento, ou declaram
explicitamente que apenas K varia com Gamma fixo.

## Fermi degenerado de nêutrons livres

Com x=p_F/(m_n c), g=2 e B=m_n⁴c⁵/(pi² hbar³):

\[
n_b=\frac{(m_nc)^3x^3}{3\pi^2\hbar^3},\quad\rho_0=m_n n_b,
\]
\[
\varepsilon=\frac B8[x\sqrt{1+x^2}(1+2x^2)-\operatorname{asinh}x],
\]
\[
P=\frac B{24}[x\sqrt{1+x^2}(2x^2-3)+3\operatorname{asinh}x].
\]

h=ln(sqrt(1+x²)); c_s²/c²=x²/[3(1+x²)]. O programa usa séries em x<0.08
para evitar cancelamento e verifica as fórmulas com integração independente
dos momentos. Nos limites não relativístico e ultrarrelativístico, recupera
P∝rho0^(5/3), epsilon≈rho0 c² e P≈epsilon/3, respectivamente.

Não inclui forças nucleares, prótons, elétrons, múons, equilíbrio beta,
superfluidez, crosta nem fases exóticas. A massa máxima deste modelo é uma
referência histórica computacional; não é uma previsão da massa máxima de
estrelas de nêutrons reais.

## Politrópica por trechos

P=K_i rho0^Gamma_i; epsilon=(1+a_i)rho0 c²+P/(Gamma_i-1).
K_{i+1}=K_i rho_break^(Gamma_i-Gamma_{i+1}).
a_{i+1}=a_i+P_break/(rho_break c²)[1/(Gamma_i-1)-1/(Gamma_{i+1}-1)].
As constantes fixam continuidade de P, epsilon e (epsilon+P)/rho0; a derivada
e a velocidade do som podem saltar. Isso não representa uma transição de
primeira ordem. É a forma de parametrização discutida por Read et al., não
uma implementação de seus conjuntos de parâmetros publicados. Todos os
exemplos internos são identificados como sintéticos e não calibrados.

## Velocidade do som constante / matéria autoligada

P=s(epsilon-epsilon_s), 0<s<=1; epsilon_s>0.
h=s/(1+s) ln[1+(1+s)P/(s epsilon_s)] e
rho0=rho0_s exp(h/s). rho0_s é uma normalização bariônica INDEPENDENTE,
obrigatória na entrada. A pressão zera com densidade de energia finita.
Não se presume uma composição de quarks nem uma fase nuclear específica.
O código não une essa EOS a matéria hadrônica através de um salto.

## EOS incompressível

epsilon=constante, h=ln(1+P/epsilon). O som tem velocidade infinita. Serve
somente para comparação com a solução interior exata de Schwarzschild.
O usuário deve desativar explicitamente a verificação de causalidade.
Não se fornece rho0 física nem massa bariônica para esse caso.

## Tabelas externas

O importador recebe uma barótropa fria de 1 dimensão; não lê diretamente uma
EOS de temperatura/composição arbitrárias. Requer pelo menos oito linhas,
P, epsilon e rho0 positivos, finitos e estritamente crescentes. Não ordena
nem elimina linhas silenciosamente. Não extrapola além da maior densidade.

Interpolação PCHIP em log epsilon(log P) e log rho0(log P) preserva
monotonicidade das quantidades, mas não impõe a primeira lei exatamente.
A primeira lei é verificada em uma malha de pontos; o erro máximo e o
limite de aceitação (5% por padrão, ajustável na API) ficam nos metadados.
5% é um filtro de importação, NÃO uma precisão científica recomendada.
Para publicação, use tabela densa, examine o resíduo e refine a malha.

A entalpia integra P/(epsilon+P) em log P com Simpson cumulativo; a inversa
log P(h) usa PCHIP. Há erro adicional de interpolação da entalpia, registrado
pela identidade dP/dh=epsilon+P. Testes variam tanto o número de linhas da
EOS quanto a malha de entalpia; reduzir rtol não elimina esses erros.

A tabela deve ter energia total, incluindo repouso. Um arquivo de energia
interna apenas pode passar na monotonicidade e ainda ser fisicamente errado;
por isso é obrigatório auditar a convenção do fornecedor. O programa não
“descobre” composição ou convenções a partir de números.

Como tabelas frias normalmente não vão até P=0, exige-se tail_gamma. Abaixo
da primeira linha aplica-se uma extensão politrópica declarada. Seu coeficiente
de repouso é ajustado para continuidade de epsilon e potencial químico.
Essa extensão pode alterar massa/raio e NÃO equivale a uma crosta nuclear
validada. Compare pontos de junção e crostas físicas antes de usar resultados
para inferência observacional. Não há extensão automática em alta densidade.

Veja data/README.md para esquemas SI, CGS e nucleares, unidades de CompOSE,
proveniência, hash SHA-256 e procedimento para obter tabelas científicas.
