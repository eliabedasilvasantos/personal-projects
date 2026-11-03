# Desempenho e expansões computacionais

A prioridade desta versão é correção e rastreabilidade. As integrações
individuais têm apenas três variáveis; o custo importante é repetir a EOS
durante muitas estrelas. NumPy organiza perfis e malhas; SciPy controla
integração adaptativa, interpolação, quadratura e otimização local.

## Implementado e testado

Varreduras com ProcessPoolExecutor e contexto spawn, preservando a ordem.
workers=1 é padrão. O teste compara sete densidades em execução serial e
com dois processos; M e R coincidem dentro do limite de 1e-14 relativo.
Processos têm custo de inicialização; duas ou três estrelas podem ficar
mais lentas em paralelo. Meça em sua máquina antes de aumentar workers.

Cada tarefa recebe EOS/configuração serializáveis. Não há cache global
mutável. Perfis são amostrados depois da integração: pedir mais amostras não
exige passos físicos menores. O passo máximo continua sendo um controle
independente. Os arrays de saída do integrador são somente leitura.

## Numba: expansão preparada, ainda não ativada

O primeiro candidato é o cálculo escalar dos integrais de Fermi e RHS
geométrico, hoje isolados das exportações e da UI. Compilar apenas uma
função dentro de um callback Python pode não reduzir significativamente o
custo total. Não compile solve_ivp ou PCHIP genericamente e não use fastmath
sem repetir os benchmarks e a checagem de causalidade.

Procedimento futuro: medir perfil de tempo → adicionar kernel opcional com
mesmas unidades/contrato → validar limites e quadratura independente →
comparar M–R, evento superficial e tolerâncias → ativar somente se houver
ganho mensurado. Numba não é dependência nem desempenho alegado nesta entrega.

## JAX: alternativa experimental futura

SciPy solve_ivp não é automaticamente diferenciável ao colocar um jax.jit
no callback. Uma implementação JAX exigirá integrador adequado, precisão
64 bits explicitamente habilitada, eventos e EOS compatíveis. Comparações
devem separar erro de integração, interpolação e diferenciação do local da
superfície. Quebras de derivadas de EOS por trechos e mudanças de estabilidade
exigem cuidado adicional. Não há backend JAX fictício no programa atual.

## Antes de HPC

Defina EOS e unidades; confirme causalidade; resolva convergência de M, R e
máximo local; refine a malha central; depois aumente paralelismo. Registre
versões, CPU, número de processos e opções de BLAS na pesquisa. Threads de
BLAS multiplicadas por processos podem causar excesso de concorrência.
Não são prometidos ganhos de velocidade sem medição.
