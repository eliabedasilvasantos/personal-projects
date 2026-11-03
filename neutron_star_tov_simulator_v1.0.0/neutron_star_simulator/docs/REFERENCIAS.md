# Referências primárias verificadas

Consultadas para esta implementação em 2 de novembro de 2026. As equações e
derivações do projeto são apresentadas em notação própria. Não há cópia
integral de artigos nem dados observacionais redistribuídos.

1. J. R. Oppenheimer; G. M. Volkoff (1939). **On Massive Neutron Cores**.
   Physical Review 55, 374–381. [DOI 10.1103/PhysRev.55.374](https://doi.org/10.1103/PhysRev.55.374).
   [Cópia do artigo original, IST Lisboa](https://blackholes.tecnico.ulisboa.pt/gritting/pdf/gravity_and_general_relativity/Oppenheimer-Volkoff_On-Massive-Neutron-Cores.pdf).
   Usado: EOS livre relativística, variável t e quatro casos finitos da Tabela I.
   Os valores históricos são arredondados; não são uma referência de alta precisão.

2. K. Schwarzschild (1916), tradução de S. Antoci (1999).
   **On the Gravitational Field of a Sphere of Incompressible Fluid According
   to Einstein's Theory**. [arXiv:physics/9912033](https://arxiv.org/abs/physics/9912033).
   Usado: solução interior exata de densidade de energia constante. Benchmark
   matemático acausal; não é microfísica nuclear.

3. J. S. Read; B. D. Lackey; B. J. Owen; J. L. Friedman (2009).
   **Constraints on a phenomenologically parametrized neutron-star equation
   of state**. Physical Review D 79, 124032.
   [DOI 10.1103/PhysRevD.79.124032](https://doi.org/10.1103/PhysRevD.79.124032),
   [arXiv:0812.2163](https://arxiv.org/abs/0812.2163).
   Usado: forma de parametrização por politrópicas contínuas. Os parâmetros
   de EOS nomeadas do artigo NÃO foram transcritos nem alegados neste projeto.

4. M. G. Alford; S. Han; M. Prakash (2013). **Generic conditions for stable
   hybrid stars**. Physical Review D 88, 083013.
   [DOI 10.1103/PhysRevD.88.083013](https://doi.org/10.1103/PhysRevD.88.083013),
   [arXiv:1302.4732](https://arxiv.org/abs/1302.4732).
   Contexto da parametrização CSS. O programa usa somente a relação linear
   como EOS autoligada de teste; não implementa a transição híbrida desse artigo.

5. J. L. Friedman; J. R. Ipser; R. D. Sorkin (1988). **Turning-point method
   for axisymmetric stability of rotating relativistic stars**.
   Astrophysical Journal 325, 722–724.
   [DOI 10.1086/166043](https://doi.org/10.1086/166043),
   [artigo original ADS](https://adsabs.harvard.edu/pdf/1988ApJ...325..722F).
   Contexto e hipóteses dos critérios de retorno; o teorema rotativo/secular
   não é aplicado indiscriminadamente como uma solução de modos dinâmicos.

6. S. Typel et al. (2022). **CompOSE reference manual**. European Physical
   Journal A 58, 221.
   [DOI 10.1140/epja/s10050-022-00847-y](https://doi.org/10.1140/epja/s10050-022-00847-y).
   [CompOSE](https://compose.obspm.fr/). Usado: convenções nucleares, condições
   de tabelas, extração e disponibilidade de M–R. Ao usar uma EOS externa,
   cite também seus autores, não apenas a base CompOSE.

7. NIST, **2022 CODATA recommended values**, publicação de dados de 2024.
   [Banco oficial](https://physics.nist.gov/cuu/Constants/),
   [listagem oficial](https://physics.nist.gov/cuu/pdf/all.pdf).
   Usado: G e m_n; c, h e carga elementar são constantes SI definidas.

8. E. E. Mamajek et al. (2015). **IAU 2015 Resolution B3 on Recommended
   Nominal Conversion Constants for Selected Solar and Planetary Properties**.
   [arXiv:1510.07674](https://arxiv.org/abs/1510.07674).
   Usado: parâmetro solar nominal (GM_sun)^N, sem alegar uma massa solar exata.

9. SciPy, documentação oficial de
   [solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html).
   Usado: métodos adaptativos, tolerâncias por componente, eventos e saída
   densa. Versão efetivamente executada registrada em cada manifest.json.

O limite newtoniano n=1 é verificado pela solução elementar theta=sin(xi)/xi,
derivada no teste; não depende de números supostamente extraídos de uma
fonte não identificada. Nenhum resultado de observação é fabricado para
“validar” os exemplos sintéticos.
