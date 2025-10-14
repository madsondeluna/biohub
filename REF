# I. Motivação
A análise de proteínas é um elemento central em bioinformática estrutural, pois permite compreender propriedades físico-químicas e estruturais que impactam diretamente na função biológica e em potenciais aplicações biotecnológicas. Atualmente, existem ferramentas que realizam análises isoladas, mas muitas vezes não são integradas, dificultando a rotina de pesquisadores que precisam lidar com dados de sequências primárias (em formato FASTA) e estruturas tridimensionais (formatos PDB e mmCIF/PDBx).

Neste projeto, o grupo propõe o desenvolvimento de um Hub para centralizar ferramentas de bioinformática em uma única plataforma, acessível via interface web, e por execução local de scripts. Essa abordagem busca reduzir a fragmentação das ferramentas, otimizar o fluxo de trabalho de análise e fornecer resultados reprodutíveis, visualmente organizados e interpretáveis.

Um diferencial importante do BioHub é sua implementação utilizando apenas a biblioteca padrão do Python 3, sem dependências externas para as funcionalidades principais. Esta escolha metodológica garante portabilidade, facilita a manutenção e torna a ferramenta acessível em qualquer ambiente computacional, desde laptops até servidores de alto desempenho, eliminando problemas comuns de incompatibilidade de versões e gerenciamento de dependências.
 
# II. Referencial Teórico
A caracterização físico-química de proteínas, a partir tanto de suas sequências primárias quanto de suas estruturas tridimensionais, é uma etapa fundamental na bioinformática e na biologia molecular. Métricas como o ponto isoelétrico (pI) permitem determinar o pH no qual a proteína apresenta carga líquida zero, informação essencial para otimizar processos de purificação, solubilidade e estabilidade em diferentes condições experimentais. A carga total e a hidrofobicidade, por sua vez, fornecem um panorama da distribuição de cargas e da organização de regiões hidrofóbicas, o que influencia diretamente a formação de domínios estruturais e a capacidade da proteína de interagir com outras moléculas ou membranas biológicas. O índice GRAVY (do inglês, Grand Average of Hydropathicity) (Kyte e Doolittle, 1982) complementa essa análise ao inferir as propriedades de solubilidade e tendência de agregação a partir do cálculo de caráter hidrofóbico ou hidrofílico global da proteína, aspecto relevante para inferir propriedades de solubilidade e tendência de agregação.

Um outro  método descrito por Kyte-Doolittle amplamente utilizado é capaz de gerar perfis de hidrofobicidade com base nos resíduos de aminoácidos para a identificação de segmentos candidatos a hélices transmembranares, assim, auxiliando na predição da topologia de proteínas de membrana (Kyte e Doolittle, 1982). Além disso, a incorporação de métricas estruturais amplia a compreensão funcional: o cálculo da superfície acessível ao solvente (SASA, em inglês) (Lee e Richards, 1971) revela regiões expostas e enterradas, fornecendo indícios sobre estabilidade conformacional e potenciais sítios de interação com ligantes. De modo complementar, a análise eletrostática por meio do APBS (do inglês Adaptive Poisson-Boltzmann Solver) (Jurrus et al., 2018) permite avaliar a distribuição de potenciais de carga na superfície da proteína, crucial para compreender mecanismos de reconhecimento molecular e afinidade de ligação.

Embora existam ferramentas consolidadas com o propósito de realizar tais análises, como o ExPASy (Gasteiger et al., 2003), ProtParam (Gasteiger et al., 2005), PyMOL (Schrödinger e DeLano, 2020), DSSP (Kabsch e Sander, 1983) e o próprio APBS (Jurrus et al., 2018), elas frequentemente operam de maneira fragmentada, exigindo que o usuário transite entre múltiplos programas e formatos de arquivo. Nesse sentido, a integração dessas métricas em uma plataforma unificada tem o potencial de facilitar a caracterização de proteínas de forma mais ágil, acessível e reprodutível, consolidando em um único ambiente análises que atualmente demandam tempo e conhecimento técnico disperso.

# III. Metodologia
### 3.1 Módulos Funcionais

As funções integradas foram representadas em oito módulos, cada um responsável por uma categoria específica de análise bioinformática. Essa implementação garante independência funcional, onde falhas em módulos individuais não comprometem a operacionalidade global do sistema.

**Módulo fetchpdb:** Este módulo automatiza a obtenção de estruturas cristalográficas diretamente do repositório RCSB PDB através de requisições HTTP. Além do download da estrutura tridimensional, o módulo executa parsing automático de informações de cabeçalho (HEADER, TITLE, COMPND, SOURCE), extraindo metadados relevantes como título da estrutura, data de deposição, classificação funcional e organismo de origem, proporcionando contexto biológico para análises subsequentes.

**Módulo fasta:** Responsável pela extração de sequências primárias a partir de coordenadas atômicas em arquivos PDB, este módulo processa sistematicamente registros ATOM, identifica a primeira cadeia polipeptídica presente na estrutura, rastreia números de resíduos para evitar duplicações, e aplica conversão de códigos de três letras para notação de letra única, gerando saída em formato FASTA padronizado compatível com ferramentas de análise de sequências.

**Módulo csv2fasta:** Este módulo oferece funcionalidade de conversão de dados tabulares para formato FASTA, facilitando a integração de dados provenientes de planilhas, bancos de dados relacionais ou resultados de experimentos de alto rendimento. O módulo suporta delimitadores customizados (vírgula, ponto-e-vírgula, tabulação), permite mapeamento flexível de colunas através de índices numéricos ou nomes de cabeçalho, e processa automaticamente múltiplas sequências em operação de lote.

**Módulo physchem:** Constitui o núcleo analítico da plataforma para caracterização físico-química de sequências proteicas. O módulo calcula peso molecular aplicando correção para perda de moléculas de água nas ligações peptídicas; determina o ponto isoelétrico através de algoritmo iterativo que resolve a equação de Henderson-Hasselbalch considerando grupos ionizáveis terminais e cadeias laterais; computa o índice GRAVY (Grand Average of Hydropathicity) como indicador do caráter hidrofóbico global; estima o coeficiente de extinção molar a 280 nm baseado em resíduos aromáticos; calcula o índice alifático correlacionado à estabilidade térmica; quantifica o índice de instabilidade baseado em frequências de dipeptídeos; estima a meia-vida celular segundo a regra N-end; e fornece composição aminoacídica detalhada com classificação em resíduos ácidos, básicos, polares e apolares.

**Módulo contacts:** Este módulo identifica interações intramoleculares através do cálculo sistemático de distâncias euclidianas entre átomos de carbono alfa (Cα). O algoritmo gera uma matriz de distâncias completa, filtra possíveis contatos baseados em um limiar configurável (usando como padrão a distância de 8.0 Å), exclui interações de resíduos sequencialmente adjacentes, e produz uma lista ordenada de pares de resíduos em contato, permitindo a identificação de núcleos hidrofóbicos, padrões de enovelamento e interações estabilizadoras de longo alcance.

**Módulo exposure:** Implementa o método de Kyte-Doolittle com janela deslizante para predição de regiões expostas ao solvente versus regiões enterradas. O algoritmo aplica uma janela móvel de tamanho configurável (usando como padrão o tamanho de 9 resíduos), calcula média de hidropaticidade para cada posição, processa bordas com janelas adaptativas, e classifica regiões em hidrofóbicas (provavelmente enterradas) ou hidrofílicas (provavelmente expostas), técnica fundamental para análise de topologia de proteínas de membrana e identificação de hélices transmembrana.

**Módulo sasa:** Representa uma implementação completa do algoritmo de Shrake-Rupley para cálculo de área de superfície acessível ao solvente. O módulo gera distribuição uniforme de pontos na esfera atômica através do método de Fibonacci (padrão 960 pontos), aplica raio de sonda configurável (padrão 1.4 Å para moléculas de água), executa teste geométrico de oclusão por átomos vizinhos para cada ponto de teste, calcula área acessível por átomo através de fração de pontos não oclusos, e agrega resultados por resíduo, fornecendo tanto SASA total da molécula quanto contribuições individuais por resíduo.

**Módulo apbs:** Funciona como wrapper automatizado para o pipeline completo de cálculo eletrostático através da equação de Poisson-Boltzmann. O módulo integra PDB2PQR para adição de hidrogênios e atribuição de cargas parciais baseadas em campo de força AMBER; gera arquivos de configuração APBS com parâmetros de grid otimizados; executa cálculos de energia de solvatação em meio aquoso e meio de referência; extrai diferenças de energia através de parsing da saída; e fornece valor final de energia de solvatação eletrostática em kJ/mol, com opção de preservação de arquivos intermediários para análises adicionais ou depuração.

# IV. Referências
GASTEIGER, Elisabeth et al. ExPASy: the proteomics server for in-depth protein knowledge and analysis. Nucleic Acids Research, v. 31, n. 13, p. 3784–3788, 1 jul. 2003.

GASTEIGER, Elisabeth et al. Protein Identification and Analysis Tools on the ExPASy Server. In: WALKER, John M. (Org.). The Proteomics Protocols Handbook. Totowa, NJ: Humana Press, 2005. p. 571–607.

JURRUS, Elizabeth et al. Improvements to the APBS biomolecular solvation software suite. Protein Science, v. 27, n. 1, p. 112–128, 2018.

KABSCH, Wolfgang; SANDER, Christian. Dictionary of protein secondary structure: Pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, v. 22, n. 12, p. 2577–2637, 1983.

KYTE, Jack; DOOLITTLE, Russell F. A simple method for displaying the hydropathic character of a protein. Journal of Molecular Biology, v. 157, n. 1, p. 105–132, 5 maio 1982.

LEE, B.; RICHARDS, F. M. The interpretation of protein structures: Estimation of static accessibility. Journal of Molecular Biology, v. 55, n. 3, p. 379-IN4, 14 fev. 1971.

SCHRÖDINGER, L.; DELANO, W. PyMOL. Disponível em: https://www.scirp.org/reference/referencespapers?referenceid=3860375. Acesso em: 14 set. 2025.

SHRAKE, A.; RUPLEY, J. A. Environment and exposure to solvent of protein atoms. Lysozyme and insulin. Journal of Molecular Biology, v. 79, n. 2, p. 351–371, 1973.
 








