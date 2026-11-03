# Importação de EOS científicas

O arquivo synthetic_polytrope_si.csv foi gerado da EOS analítica
P=0.01 rho0², epsilon=rho0 c²+P, Gamma=2. Não contém dados experimentais,
observacionais nem uma EOS nuclear publicada. É somente uma fixture e um
exemplo de formato executável.

## Formatos canônicos aceitos

CSV com cabeçalho, vírgula como separador, ponto decimal e notação científica.
Pelo menos oito linhas. Não inclua comentários, unidades misturadas, linhas
duplicadas, NaN ou zero. Linhas devem ter densidade, pressão e energia total
estritamente crescentes. A extensão até P=0 será declarada separadamente.

| --table-units | Colunas exatas | Convenção |
|---|---|---|
| si | rho0_kg_m3,pressure_pa,epsilon_j_m3 | kg/m³, Pa, J/m³ |
| cgs | rho0_g_cm3,pressure_dyn_cm2,epsilon_erg_cm3 | g/cm³, dyn/cm², erg/cm³ |
| nuclear | n_b_fm3,pressure_mev_fm3,epsilon_mev_fm3 | fm^-3, MeV/fm³, MeV/fm³ |

No formato nuclear, --baryon-mass-kg é obrigatório para rho0=m_b n_b.
A energia epsilon permanece a energia TOTAL fornecida pela EOS; m_b não é
usado para substituir/restaurar um termo de repouso que já está nessa energia.
Para CGS: 1 erg/cm³=0.1 J/m³. Não confunda energia por bárion (MeV) com
energia por volume (MeV/fm³).

## Como obter uma fonte científica

1. Acesse [CompOSE](https://compose.obspm.fr/) e leia a ficha do modelo
   escolhido e o [manual oficial](https://doi.org/10.1140/epja/s10050-022-00847-y).
   Selecione matéria de estrela de nêutrons fria em equilíbrio beta e
   neutralidade de carga quando essas forem as hipóteses desejadas. “Neutron
   matter” puro e matéria simétrica não são equivalentes a esse caso.
2. Confira se o modelo tem T=0 ou somente temperatura mínima finita, quais
   espécies/leptons estão incluídos, domínio em densidade, definição de energia,
   convenção bariônica e descrição da crosta. Uma fatia de temperatura fixa
   de uma tabela geral não se torna uma EOS fria por renomear o arquivo.
3. Baixe os arquivos e a ficha, registre versão/data, publicações dos autores
   e licença. Use o software oficial CompOSE para extrair uma sequência 1D
   nas condições escolhidas e as quantidades n_b, P e epsilon total. O pacote
   não inventa uma conversão universal das colunas cruas de eos.thermo.
4. Exporte para um dos CSVs acima. Se a fonte fornece energia total por bárion
   e_tot, epsilon=n_b e_tot nas unidades nucleares correspondentes. Se fornece
   apenas energia interna por bárion, é necessário aplicar a convenção de
   repouso exata da fonte; não some m_n c² indiscriminadamente. Audite pelo
   menos uma linha manualmente com os dados da ficha.
5. Informe --provenance com modelo, DOI, versão, composição, T e convenções.
   Informe --tail-gamma conscientemente: essa cauda é uma aproximação declarada,
   não uma crosta nuclear validada. EOS com superfície finita autoligada devem
   usar um tratamento próprio; a cauda deste importador é de densidade nula.
6. Examine o resíduo da primeira lei, a derivada da entalpia, c_s²/c², ausência
   de extrapolação, sensibilidade ao ponto de junção/cauda, número de linhas
   e tolerâncias. Compare M–R com o arquivo eos.mr do fornecedor, **quando
   disponível**, usando as mesmas condições e crosta. Registre discrepâncias.

## O que é recusado

Saltos de densidade, trechos de pressão constante (transição de primeira
ordem), comprimentos incompatíveis, energia/pressão não monotônicas e domínio
insuficiente em alta densidade. O importador não suaviza saltos de fase nem
elimina linhas para fazer a integração passar. EOS autoligadas tabeladas não
são suportadas pela cauda atual. Também não reconstrói uma fatia em equilíbrio
beta a partir de uma tabela multidimensional.

Os CSVs externos não são redistribuídos nesta entrega. O hash do arquivo
efetivamente usado fica nos metadados da EOS e no manifesto de execução.
O checksum garante identidade dos bytes; não garante validade científica.
