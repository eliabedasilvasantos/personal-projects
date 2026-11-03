# Fundamentação física e numérica

## 1. Problema físico e hipóteses

Equilíbrio estático, esférico, assintoticamente plano de um fluido perfeito
isotrópico em relatividade geral, sem rotação, campo magnético, anisotropia,
elasticidade, transporte térmico ou evolução temporal. A EOS é uma barótropa
fria, com composição já incorporada. Resolver TOV não calcula a composição
nuclear nem demonstra estabilidade dinâmica.

Usa-se o raio areal r e a métrica

\[
ds^2=-e^{2\Phi(r)}c^2dt^2+[1-2Gm(r)/(rc^2)]^{-1}dr^2+r^2d\Omega^2.
\]

Da equação de Einstein e da conservação covariante do tensor energia-momento:

\[
\frac{dm}{dr}=4\pi r^2\frac{\varepsilon}{c^2},\quad
\frac{dP}{dr}=-\frac{G(\varepsilon+P)(m+4\pi r^3P/c^2)}
{c^2r^2[1-2Gm/(rc^2)]}.
\]

Aqui epsilon é **energia total por volume próprio**, incluindo repouso e
energia interna (J/m³), P está em Pa=J/m³ e m em kg. A densidade que entra em
dm/dr é epsilon/c², não a densidade de repouso rho0. Esta última é m_b n_b e
depende da convenção explícita para a massa bariônica m_b.

Dimensões: dm/dr tem kg/m e dP/dr tem kg/(m² s²). Quando P << epsilon,
2Gm/(rc²) << 1 e 4pi r³P/c² << m, recupera-se dP/dr=-G rho m/r².
O código integra as equações relativísticas completas.

## 2. Unidades e adimensionalização exata

Para uma escala L, x=r/L, mu=Gm/(c²L), e=G L² epsilon/c⁴,
p=G L² P/c⁴ e b=G L² rho0/c². Então

\[
\mu'=4\pi x^2 e,\quad
p'=-(e+p)(\mu+4\pi x^3p)/(x^2(1-2\mu/x)).
\]

Não se elimina c de uma grandeza SI sem converter sua dimensão. A escala L
é numérica, não um parâmetro físico: variar L deve preservar M e R.
1 MeV/fm³ = 1.602176634e32 J/m³ = 1.602176634e32 Pa;
1 g/cm³ = 1000 kg/m³; 1 dyn/cm²=0.1 Pa; 1 fm^-3=1e45 m^-3.
Em unidades nucleares hbar=c=1, é ainda necessário declarar se a entrada é
energia por bárion, energia por volume ou massa por volume. O importador
aceita energia TOTAL por volume, jamais energia interna não identificada.

G e m_n seguem CODATA 2022. c, h e a carga elementar são definidos no SI.
Usa-se (GM_sun)^N=1.3271244e20 m³/s² (IAU 2015 B3); a massa solar em kg é
definida aqui por esse parâmetro nominal dividido por G. Não se propagam
incertezas dessas constantes nesta versão.

## 3. Entalpia e superfície

\[
h(P)=\int_0^P\frac{dP'}{\varepsilon(P')+P'},\quad
\frac{dh}{dx}=-\frac{\mu+4\pi x^3p}{x^2(1-2\mu/x)}.
\]

Integram-se (mu,h,mu_b) com solve_ivp, raio como variável independente.
O evento h=0, com direção negativa, localiza P=0. Evita-se uma pressão
mínima artificial e o mau condicionamento da raiz em P quando P' tende a
zero. A EOS só é avaliada em max(h,0) nas etapas provisórias do método;
esses valores exteriores são extensões numéricas e não entram no perfil
físico, que termina na raiz. Em matéria autoligada a densidade superficial
finita corresponde ao limite **interior**; o exterior é vácuo.

As tolerâncias absolutas são para mu, h e mu_b, adimensionais, e não para
kg/Pa. DOP853 é a opção padrão; RK45/Radau/BDF/LSODA também estão expostos.
Refinar tolerâncias, passo máximo, escala L e raio inicial é obrigatório
em qualquer novo regime. A raiz é interpolada pelo integrador, não pela
última amostra de um gráfico. Falha do integrador, horizonte, causalidade
violada e superfície fora do domínio produzem exceções explícitas.

## 4. Centro regular

Se e_c e p_c são finitos, com A=2pi(e_c/3+p_c),

\[
\mu(x)=4\pi e_c x^3/3+O(x^5),\quad
h(x)=h_c-Ax^2+O(x^4),\quad
\mu_b(x)=4\pi b_c x^3/3+O(x^5).
\]

O integrador começa em r0>0 com essa série. Exige-se correção central
pequena diante de h_c. A amostra r=0 é preenchida pelas condições exatas.
Não se substitui r por um epsilon arbitrário no denominador da TOV.

## 5. Massa gravitacional, bariônica e observáveis

M=m(R) é a massa gravitacional total do espaço-tempo estático. Não há uma
segunda “massa total” independente. A massa bariônica (se uma convenção de
rho0 é conhecida) é

\[
M_b=\int_0^R4\pi r^2\rho_0(1-2Gm/(rc^2))^{-1/2}dr.
\]

Ela é uma contagem de repouso em volume próprio e difere da integral de
epsilon em coordenada areal. Energia de ligação referida a essa convenção:
(M_b-M)c². A compactação relatada é C=GM/(Rc²); o fator de Schwarzschild é
1-2C. Redshift de um fóton superficial a um observador estático no infinito:
z=(1-2C)^(-1/2)-1. No centro, Gm/(rc²) tem limite zero.

## 6. Estabilidade: alcance da inferência

Ao longo de uma sequência suave com a MESMA EOS fria e sem rotação, um
primeiro máximo de M(epsilon_c), conectado a um ramo de baixa densidade
radialmente estável, sinaliza o início da instabilidade radial sob as
hipóteses usuais. dM/depsilon_c é um indicador condicionado a essas
hipóteses. Um ramo com derivada positiva após extremos adicionais não é
automaticamente estável: modos instáveis anteriores precisam ser contados.

As classificações do programa são “candidato no ramo inicial”, “derivada
negativa”, “inconclusivo após extremos”, e “região de extremo”. O programa
não resolve o problema de autovalores de pulsações radiais omega². Não
equaciona equilíbrio hidrostático com estabilidade dinâmica. Transições de
primeira ordem com saltos de energia não são suportadas pelo importador;
exigem tratamento de interfaces e análise de estabilidade específica.

Um máximo refinado é local e delimitado pela malha explorada. Se o maior
valor estiver na borda, não se anuncia uma massa máxima global. Limites
do intervalo e erros de modelos não são incertezas estatísticas.
