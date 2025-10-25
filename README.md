# BioHub: Uma Plataforma para Análise de Sequências e Estruturas de Proteínas

![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Research-006699?logo=google-scholar&logoColor=white)
[![Python](https://img.shields.io/badge/Python-3.x-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![CLI](https://img.shields.io/badge/Command%20Line-CLI-orange?logo=gnometerminal&logoColor=white)](#)
[![PDB](https://img.shields.io/badge/Protein%20Data%20Bank-PDB-0052CC?logo=databricks&logoColor=white)](https://www.wwpdb.org/)
[![CSV](https://img.shields.io/badge/Data-CSV-FFD700?logo=files&logoColor=black)](https://en.wikipedia.org/wiki/Comma-separated_values)
[![FASTA](https://img.shields.io/badge/Sequence-FASTA-808080?logo=dna&logoColor=white)](https://en.wikipedia.org/wiki/FASTA_format)




**BioHub** é uma ferramenta de linha de comando leve, prática e centralizadora, escrita em Python puro, para realizar análises bioinformáticas a partir de arquivos de estrutura de proteínas (PDB) e sequências de aminoácidos.

---

```
 _._     _,-'""`-._
(,-.`._,'(       |\`-/|
    `-.-' \ )-`( , o o)
          `-    \`_`"'-
█████████████████████████████████████████████████████████████████████
█▌                                                                 ▐█
█▌                                                                 ▐█
█▌    .______    __    ______    __    __   __    __  .______      ▐█
█▌    |   _  \  |  |  /  __  \  |  |  |  | |  |  |  | |   _  \     ▐█
█▌    |  |_)  | |  | |  |  |  | |  |__|  | |  |  |  | |  |_)  |    ▐█
█▌    |   _  <  |  | |  |  |  | |   __   | |  |  |  | |   _  <     ▐█
█▌    |  |_)  | |  | |  `--'  | |  |  |  | |  `--'  | |  |_)  |    ▐█
█▌    |______/  |__|  \______/  |__|  |__|  \______/  |______/     ▐█
█▌                                                                 ▐█
█▌                                                                 ▐█
█████████████████████████████████████████████████████████████████████
```

---

## Equipe Desenvolvedora (Programa de Pós Graduação em Bioinformática - UFMG)

- Alisson Clementino da Silva 
- Ana Karoline da Nóbrega Nunes Alves
- Laura da Silva Ribeiro de Souza
- Leonardo Henrique da Silva
- Madson Allan de Luna Aragão

---

## Visão Geral

Este projeto foi desenvolvido com o objetivo de fornecer uma solução simples, rápida e altamente portátil para tarefas comuns na análise de proteínas. É ideal para:

* **Fins educacionais**: Demonstra a implementação de algoritmos bioinformáticos fundamentais sem a abstração de bibliotecas complexas.
* **Integração em pipelines**: A natureza leve e a interface de linha de comando facilitam a automação de análises em scripts de shell.
* **Acessibilidade**: Por não requerer bibliotecas externas, a ferramenta pode ser executada em qualquer ambiente com uma instalação padrão do Python 3, eliminando problemas de gerenciamento de dependências.

---

## Como Funciona

A BioHub foi intencionalmente construída utilizando apenas a biblioteca padrão do Python. Constantes físico-químicas, como pesos moleculares e escalas de hidropatia, são armazenadas localmente em dicionários.

**Justificativa**: Esta abordagem garante a portabilidade e simplicidade para análise de estrutura de proteínas. Embora o uso de bancos de dados externos ou pacotes especializados pudesse oferecer maior precisão ou dados mais atualizados, a implementação atual é suficiente para fins demonstrativos e educacionais, focando na lógica algorítmica.

---

## Funcionalidades Detalhadas

A ferramenta é organizada em oito subcomandos principais:

### 1. `fetchpdb`

Faz o download automático de arquivos PDB diretamente do banco de dados RCSB (Protein Data Bank). Basta fornecer o código PDB de 4 caracteres e a ferramenta baixará a estrutura correspondente.

**Recursos adicionais:**
* Extrai informações do cabeçalho do arquivo PDB (título, data de deposição, classificação, organismo de origem)
* Permite especificar um nome customizado para o arquivo de saída
* Validação automática do código PDB fornecido

### 2. `fasta`

Converte um arquivo PDB em uma sequência no formato FASTA. O processo foca nos registros `ATOM` e extrai a sequência de aminoácidos da primeira cadeia (chain) encontrada no arquivo.

### 3. `csv2fasta`

Converte arquivos CSV contendo sequências de aminoácidos para o formato FASTA. Esta funcionalidade é especialmente útil para processar grandes volumes de sequências armazenadas em planilhas ou bancos de dados tabulares. O CSV input deve conter apenas duas colunas com as variáveis "ID" e "sequência de aminoácidos".

**Recursos configuráveis:**
* Suporte para delimitadores customizados (vírgula, ponto-e-vírgula, tabulação, etc.)
* Especificação de colunas por nome ou índice numérico
* Opção para processar cabeçalhos na primeira linha
* Flexibilidade para mapear qualquer coluna como identificador ou sequência

### 4. `physchem`

Calcula um conjunto expandido de propriedades físico-químicas essenciais para uma dada sequência de aminoácidos:

* **Peso Molecular (MW)**: A massa da molécula, fundamental para experimentos de espectrometria de massa e SDS-PAGE.
* **Ponto Isoelétrico (pI)**: O pH no qual a proteína tem carga líquida zero, influenciando sua solubilidade e comportamento em cromatografia de troca iônica.
* **Coeficiente de Extinção**: Estima como a proteína absorve luz a 280 nm, útil para determinar a concentração da proteína via espectrofotometria.
* **GRAVY (Grand Average of Hydropathicity)**: Um índice do caráter hidrofóbico ou hidrofílico geral de uma proteína. Valores positivos indicam hidrofobicidade (ex.: proteínas de membrana), enquanto valores negativos indicam hidrofilicidade (ex.: proteínas citosólicas).
* **Índice Alifático**: Medida do volume ocupado por cadeias laterais alifáticas (Ala, Val, Ile, Leu), correlacionado com a estabilidade térmica da proteína.
* **Índice de Instabilidade**: Prevê se a proteína é estável ou instável em ambientes fisiológicos (valores > 40 indicam instabilidade).
* **Meia-vida**: Estimativa da meia-vida da proteína em mamíferos (in vitro) com base no resíduo N-terminal.
* **Composição de Aminoácidos**: A contagem e a frequência de cada resíduo, incluindo contagens de resíduos ácidos, básicos, polares e apolares.

### 5. `contacts`

Identifica e lista contatos intramoleculares com base na distância entre os átomos de **Carbono Alfa (CA)**. A distância entre os CAs é um excelente proxy para a proximidade entre os resíduos. Este cálculo é útil para:

* Identificar o núcleo hidrofóbico (core) da proteína.
* Analisar a topologia do enovelamento.
* Estudar interações de longo alcance que estabilizam a estrutura terciária.
* Ajuste customizável do limiar de distância (padrão: 8.0 Å).

### 6. `exposure`

Calcula a **hidrofobicidade** de cada átomo da proteína utilizando a **escala de Kyte-Doolittle**. Todos os átomos de um resíduo recebem o valor de hidrofobicidade característico daquele aminoácido. Esta análise é fundamental para:

* Identificar regiões hidrofóbicas (núcleo da proteína) vs. hidrofílicas (superfície)
* Predizer exposição ao solvente baseada no caráter químico dos resíduos
* Mapear o perfil de hidrofobicidade ao longo da estrutura 3D
* Visualizar propriedades físico-químicas diretamente na estrutura (usando B-factor)

### 7. `sasa`

Calcula a **Área de Superfície Acessível ao Solvente (SASA)** usando o algoritmo de Shrake-Rupley implementado em Python puro. Esta métrica é fundamental para estudos de enovelamento, estabilidade e interações moleculares.

**Características do cálculo:**
* Raio da sonda customizável (padrão: 1.4 Å para água)
* Número de pontos ajustável para controle de precisão (padrão: 960 pontos)
* Cálculo por resíduo e SASA total da molécula
* Uso de raios de Van der Waals específicos para cada tipo de átomo

### 8. `apbs` [BETA, ainda em validação]

Calcula a **energia de solvatação eletrostática** através da resolução da equação de Poisson-Boltzmann. Este comando atua como um wrapper que automatiza todo o fluxo de trabalho:

**Pipeline automatizado:**
1. Converte o PDB para formato PQR usando PDB2PQR (adição de hidrogênios e atribuição de cargas)
2. Gera arquivo de configuração para APBS
3. Executa o cálculo eletrostático
4. Extrai e exibe a energia de solvatação final

**Requisitos:**
* PDB2PQR instalado e disponível no PATH
* APBS instalado e disponível no PATH

**Opções:**
* `--no-cleanup`: Mantém arquivos intermediários (PQR, apbs.in, etc.) para depuração ou análises adicionais

---

<p align="center">
  <img src="/img/9FC335A9-F378-42B4-80D8-35A8A3F61FF1_1_201_a.jpeg" alt="Workflow" width="100%"/>
</p>

---

## Requisitos

* Python 3.x (biblioteca padrão)

1. Verifique a instalação:

```
python3 --version
```
```
python --version
```

2. Instale python3 caso o comando não retorne a versão.

Com exceção do comando apbs que exige as dependências pdb2pqr e apbs, nenhuma outra biblioteca é necessária. 

---

## Instalação

1. Clone este repositório:

   ```bash
   git clone https://github.com/madsondeluna/biohub.git
   cd biohub
   ```

2. Torne o script executável (opcional, para conveniência):

   ```bash
   chmod +x biohub.py
   ```

---

## Uso

A ferramenta é executada a partir do terminal, seguindo o padrão:

```bash
python3 biohub.py [COMANDO] [ARGUMENTOS] [OPÇÕES]
```

Ou, se o arquivo for executável:

```bash
./biohub.py [COMANDO] [ARGUMENTOS] [OPÇÕES]
```

---

> Exemplo de tela do BioHub em uso via CLI:

<p align="center">
  <img src="/img/example.jpeg" alt="Workflow" width="100%"/>
</p>

---

### Ajuda (ou `biohub.py -h`)

```bash
python3 biohub.py -h
```

**Saída:**
```
biohub.py [-h] {fetchpdb,fasta,csv2fasta,physchem,contacts,exposure,sasa,apbs} ...

BioHub: Uma ferramenta CLI para análise de proteínas.

positional arguments:
  {fetchpdb,fasta,csv2fasta,physchem,contacts,exposure,sasa,apbs}
                        Função a ser executada
    fetchpdb            Baixa um arquivo PDB do RCSB.
    fasta               Converte um arquivo PDB em uma sequência FASTA.
    csv2fasta           Converte um arquivo CSV em um formato FASTA.
    physchem            Calcula propriedades físico-químicas de uma sequência.
    contacts            Calcula contatos intramoleculares a partir de um arquivo PDB.
    exposure            Calcula a hidrofobicidade usando a escala Kyte-Doolittle.
    sasa                Calcula a Área de Superfície Acessível ao Solvente (SASA).
    apbs[BETA]          Calcula a energia de solvatação eletrostática.

optional arguments:
  -h, --help            show this help message and exit

Use biohub.py COMANDO -h para ajuda detalhada sobre um comando.
```

---

## Exemplos Detalhados

### 1) Baixar Arquivo PDB do RCSB

Baixar a estrutura 1A2B com o nome padrão:

```bash
python3 biohub.py fetchpdb 1A2B
```

Baixar e salvar com nome customizado:

```bash
python3 biohub.py fetchpdb 1A2B -o minha_proteina.pdb
```

**Baixar apenas cadeias específicas:**

```bash
python3 biohub.py fetchpdb 4HHB --chains A,B -o hemoglobina_AB.pdb
```

**Baixar apenas proteína (remove água e ligantes):**

```bash
python3 biohub.py fetchpdb 1CRN --protein-only -o crambin_clean.pdb
```

**Combinar filtros (cadeia específica sem água/ligantes):**

```bash
python3 biohub.py fetchpdb 4HHB --chains A --protein-only -o hemo_A_clean.pdb
```

**Ajuda do comando:**
```bash
python3 biohub.py fetchpdb -h
```

**Opções disponíveis:**
* `PDB_ID`: O código de 4 caracteres do PDB (obrigatório)
* `-o, --output ARQUIVO`: Nome do arquivo de saída (padrão: `PDBID.pdb`)
* `--chains CHAINS`: Cadeias a serem mantidas, separadas por vírgula (ex: A,B). Se omitido, mantém todas
* `--protein-only`: Mantém apenas átomos de proteína (remove água, ligantes e heteroátomos)

**Benefícios da filtragem:**
- **Redução de tamanho**: Arquivos menores e mais rápidos para processar
- **Foco na análise**: Remove elementos não essenciais para análise de estrutura de proteínas
- **Preparação para simulações**: Muitos softwares de simulação requerem apenas a proteína
- **Análise por subunidades**: Permite estudar cadeias individuais de complexos proteicos

---

### 2) Converter PDB para FASTA

Exibir a sequência no terminal:

```bash
python3 biohub.py fasta proteina.pdb
```

Salvar a sequência em um arquivo:

```bash
python3 biohub.py fasta proteina.pdb -o proteina.fasta
```

**Ajuda do comando:**
```bash
python3 biohub.py fasta -h
```

**Opções disponíveis:**
* `ARQUIVO_PDB`: Caminho para o arquivo PDB de entrada (obrigatório)
* `-o, --output ARQUIVO`: Salva a saída em um arquivo FASTA (padrão: stdout)

---

### 3) Converter CSV para FASTA

Converter com configurações padrão (colunas 0 e 1, sem cabeçalho):

```bash
python3 biohub.py csv2fasta sequencias.csv
```

Converter especificando colunas por nome (com cabeçalho):

```bash
python3 biohub.py csv2fasta sequencias.csv --header --id-col "ID" --seq-col "Sequencia" -o saida.fasta
```

Converter com delimitador customizado (ponto-e-vírgula):

```bash
python3 biohub.py csv2fasta dados.csv --delimiter ";" --header --id-col 0 --seq-col 2
```

**Ajuda do comando:**
```bash
python3 biohub.py csv2fasta -h
```

**Opções disponíveis:**
* `ARQUIVO_CSV`: Caminho para o arquivo CSV de entrada (obrigatório)
* `-o, --output ARQUIVO`: Salva a saída em um arquivo FASTA (padrão: stdout)
* `--id-col COLUNA`: Coluna do identificador - nome ou índice baseado em 0 (padrão: 0)
* `--seq-col COLUNA`: Coluna da sequência - nome ou índice baseado em 0 (padrão: 1)
* `--header`: Flag para indicar que a primeira linha é um cabeçalho
* `--delimiter CHAR`: Caractere delimitador (padrão: `,`)

---

### 4) Calcular Propriedades Físico-Químicas

Analisar uma sequência fornecida diretamente:

```bash
python3 biohub.py physchem "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
```

Salvar resultados em CSV:

```bash
python3 biohub.py physchem "MKTAYIAKQ" -o propriedades.csv
```

**Ajuda do comando:**
```bash
python3 biohub.py physchem -h
```

**Opções disponíveis:**
* `SEQUENCIA`: A sequência de aminoácidos a ser analisada (obrigatório)
* `-o, --output ARQUIVO_CSV`: Salva os resultados em um arquivo CSV

**Propriedades calculadas:**
* Comprimento da sequência
* Peso molecular (Da)
* Ponto isoelétrico (pI)
* GRAVY (hidropaticidade)
* Índice alifático
* Índice de instabilidade
* Meia-vida estimada (mamíferos, in vitro)
* Total de resíduos ácidos (Asp+Glu)
* Total de resíduos básicos (Arg+Lys+His)
* Total de resíduos polares
* Total de resíduos apolares
* Composição completa de aminoácidos (contagem e percentual)

---

### 5) Calcular Contatos Intramoleculares

Encontrar contatos com o limiar padrão (8.0 Å):

```bash
python3 biohub.py contacts proteina.pdb
```

Ajustar o limiar de distância para 10.5 Å:

```bash
python3 biohub.py contacts proteina.pdb -t 10.5
```

Salvar resultados em CSV:

```bash
python3 biohub.py contacts proteina.pdb -t 8.0 -o contatos.csv
```

**Ajuda do comando:**
```bash
python3 biohub.py contacts -h
```

**Opções disponíveis:**
* `ARQUIVO_PDB`: Caminho para o arquivo PDB de entrada (obrigatório)
* `-t, --threshold FLOAT`: Distância máxima em Angstroms para considerar um contato (padrão: 8.0)
* `-o, --output ARQUIVO_CSV`: Salva os resultados em um arquivo CSV

---

### 6) Calcular Hidrofobicidade (ou Exposure)

Executar análise de hidrofobicidade por átomo:

```bash
python3 biohub.py exposure proteina.pdb
```

Salvar resultados em CSV (com detalhes por átomo):

```bash
python3 biohub.py exposure proteina.pdb -o exposicao.csv
```

**Gerar arquivo PDB anotado com valores no B-factor:**

```bash
python3 biohub.py exposure proteina.pdb --write-pdb proteina_hydro.pdb
```

**Combinar CSV e PDB anotado:**

```bash
python3 biohub.py exposure proteina.pdb -o exposicao.csv --write-pdb proteina_hydro.pdb
```

**Gerar visualização PyMOL:**

```bash
python3 biohub.py exposure proteina.pdb --write-pdb proteina_hydro.pdb --pymol proteina_hydro.pse
```

Este comando gera:
- `proteina_hydro.pdb`: PDB com hidrofobicidade no B-factor
- `proteina_hydro.pml`: Script PyMOL com visualização pré-configurada
- `proteina_hydro.pse`: Sessão PyMOL (se PyMOL estiver instalado)

**Ajuda do comando:**
```bash
python3 biohub.py exposure -h
```

**Opções disponíveis:**
* `ARQUIVO_PDB`: Caminho para o arquivo PDB de entrada (obrigatório)
* `-o, --output ARQUIVO_CSV`: Salva os resultados por átomo em CSV (Chain, ResNum, ResName, AtomNum, AtomName, Hydrophobicity)
* `--write-pdb ARQUIVO_PDB`: Gera um arquivo PDB com a hidrofobicidade escrita no B-factor de cada átomo
* `--pymol ARQUIVO_PSE`: Gera script PyMOL (.pml) e sessão (.pse) para visualização interativa

**Interpretação dos resultados:**
* **Valores positivos**: Aminoácidos hidrofóbicos (Ile, Val, Leu, Phe, etc.) - tendem a estar enterrados no núcleo da proteína
* **Valores negativos**: Aminoácidos hidrofílicos (Arg, Lys, Asp, Glu, etc.) - tendem a estar expostos ao solvente na superfície
* **Range de valores**: -4.5 (Arg, mais hidrofílico) a +4.5 (Ile, mais hidrofóbico)
* **Nota**: Todos os átomos de um mesmo resíduo recebem o mesmo valor de hidrofobicidade

**Formato de saída CSV:**
```
Chain,ResNum,ResName,AtomNum,AtomName,Hydrophobicity
A,1,THR,1,N,-0.700
A,1,THR,2,CA,-0.700
A,1,THR,3,C,-0.700
```

**Visualização do PDB anotado:**
O arquivo PDB gerado com `--write-pdb` pode ser aberto em visualizadores moleculares como PyMOL ou Chimera, onde você pode colorir a estrutura pelos valores do B-factor para visualizar o perfil de hidrofobicidade:

```bash
# No PyMOL:
load proteina_hydro.pdb
spectrum b, red_white_blue, minimum=-4.5, maximum=4.5
```

---

### 7) Calcular SASA (Área de Superfície Acessível ao Solvente)

Executar com parâmetros padrão (sonda de 1.4 Å, 960 pontos):

```bash
python3 biohub.py sasa proteina.pdb
```

Ajustar a precisão e o raio da sonda:

```bash
python3 biohub.py sasa proteina.pdb --probe-radius 1.5 --num-points 2000
```

Salvar resultados por átomo em CSV:

```bash
python3 biohub.py sasa proteina.pdb -o sasa_resultados.csv
```


**Gerar arquivo PDB anotado com SASA (média por resíduo) no B-factor:**

```bash
python3 biohub.py sasa proteina.pdb --write-pdb proteina_sasa.pdb
```

**Combinar CSV e PDB anotado:**

```bash
python3 biohub.py sasa proteina.pdb -o sasa_resultados.csv --write-pdb proteina_sasa.pdb
```

**Gerar visualização PyMOL (média por resíduo, gradiente invertido):**

```bash
python3 biohub.py sasa proteina.pdb --write-pdb proteina_sasa.pdb --pymol proteina_sasa.pse
```

Este comando gera:
- `proteina_sasa.pdb`: PDB com SASA médio por resíduo no B-factor
- `proteina_sasa.pml`: Script PyMOL com visualização pré-configurada
- `proteina_sasa.pse`: Sessão PyMOL (se PyMOL estiver instalado)

**Ajuda do comando:**

```bash
python3 biohub.py sasa -h
```

**Opções disponíveis:**

- `ARQUIVO_PDB`: Caminho para o arquivo PDB de entrada (obrigatório)
- `--probe-radius FLOAT`: Raio da sonda do solvente em Angstroms (padrão: 1.4 para água)
- `--num-points INT`: Número de pontos na superfície de cada átomo (padrão: 960)
- `-o, --output ARQUIVO_CSV`: Salva os resultados por átomo em CSV (Chain, ResNum, ResName, AtomNum, AtomName, SASA_A2)
- `--write-pdb ARQUIVO_PDB`: Gera um arquivo PDB com o SASA médio por resíduo escrito no B-factor
- `--pymol ARQUIVO_PSE`: Gera script PyMOL (.pml) e sessão (.pse) para visualização interativa

**Novidades e Observações:**

- O valor do B-factor para SASA agora é a **média por resíduo** (mais relevante biologicamente)
- O gradiente de cores para SASA foi **invertido**: vermelho = enterrado, azul = exposto
- O range de visualização é ajustado automaticamente (percentil 70 dos resíduos expostos) para maior sensibilidade
- A SASA total da molécula é sempre exibida no terminal
- Resíduos completamente enterrados terão SASA ≈ 0.00 Ų
- Resíduos totalmente expostos podem ter SASA > 10 Ų (ajustado pelo percentil)

**Formato de saída CSV:**

```csv
Chain,ResNum,ResName,AtomNum,AtomName,SASA_A2
A,1,THR,1,N,20.78
A,1,THR,2,CA,10.87
A,1,THR,3,C,0.00
```

**Visualização do PDB anotado:**

O arquivo PDB gerado com `--write-pdb` permite visualizar a acessibilidade ao solvente diretamente na estrutura 3D:

```bash
# No PyMOL:
load proteina_sasa.pdb
spectrum b, red_white_blue, minimum=0, maximum=10.8
# (range ajustado automaticamente)
```

**Exemplo prático:**

```bash
# Gerar SASA por resíduo e visualização
python3 biohub.py sasa 4HHB.pdb --write-pdb 4HHB_sasa.pdb --pymol 4HHB_sasa.pse --num-points 300

# Abrir no PyMOL
pymol 4HHB_sasa.pml
```

No PyMOL, resíduos enterrados aparecerão em vermelho, parcialmente expostos em branco, e totalmente expostos em azul.

---

### 8) Calcular Energia de Solvatação (APBS) [BETA, ainda em validação]

Executar a análise padrão (arquivos temporários são removidos):

```bash
python3 biohub.py apbs proteina.pdb
```

Executar e manter os arquivos intermediários:

```bash
python3 biohub.py apbs proteina.pdb --no-cleanup
```

**Ajuda do comando:**
```bash
python3 biohub.py apbs -h
```

**Opções disponíveis:**
* `ARQUIVO_PDB`: Caminho para o arquivo PDB de entrada (obrigatório)
* `--no-cleanup`: Previne a remoção dos arquivos temporários (PQR, apbs.in, etc.)

**Requisitos:**
* Certifique-se de que `pdb2pqr` e `apbs` estão instalados e disponíveis no PATH do sistema

**Pipeline executado:**
1. Conversão PDB → PQR (adição de hidrogênios e cargas)
2. Geração do arquivo de configuração APBS
3. Cálculo da energia de solvatação eletrostática
4. Extração e exibição do resultado final em kJ/mol

---

## Visualização com PyMOL

O BioHub pode gerar automaticamente arquivos de sessão PyMOL (.pse) e scripts (.pml) para visualização interativa dos resultados de hidrofobicidade e SASA. Esta integração permite visualizar as propriedades calculadas diretamente na estrutura 3D da proteína.

### Instalação do PyMOL

Para usar esta funcionalidade, você precisa ter o PyMOL instalado:

**macOS:**
```bash
brew install pymol
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install pymol
```

**Windows ou outras plataformas:**
Baixe de [https://pymol.org/](https://pymol.org/)

### Gerando Visualizações

#### Hidrofobicidade:

```bash
# Gera PDB anotado + script PyMOL + sessão PyMOL
python3 biohub.py exposure 1CRN.pdb --write-pdb 1CRN_hydro.pdb --pymol 1CRN_hydro.pse
```

**Arquivos gerados:**
- `1CRN_hydro.pdb`: PDB com valores de hidrofobicidade na coluna B-factor
- `1CRN_hydro.pml`: Script PyMOL com comandos de visualização
- `1CRN_hydro.pse`: Sessão PyMOL completa (se PyMOL estiver no PATH)

**Esquema de cores:**
- **Azul**: Regiões hidrofílicas (valores negativos: -4.5 a 0)
- **Branco**: Regiões neutras (próximo de 0)
- **Vermelho**: Regiões hidrofóbicas (valores positivos: 0 a +4.5)

#### SASA (Acessibilidade ao Solvente):

```bash
# Gera PDB anotado + script PyMOL + sessão PyMOL
python3 biohub.py sasa 1CRN.pdb --write-pdb 1CRN_sasa.pdb --pymol 1CRN_sasa.pse --num-points 960
```

**Arquivos gerados:**
- `1CRN_sasa.pdb`: PDB com valores de SASA na coluna B-factor
- `1CRN_sasa.pml`: Script PyMOL com comandos de visualização
- `1CRN_sasa.pse`: Sessão PyMOL completa (se PyMOL estiver no PATH)

**Esquema de cores:**
- **Azul**: Átomos enterrados (SASA baixo: 0 Ų)
- **Branco**: Parcialmente exposto
- **Vermelho**: Átomos expostos ao solvente (SASA alto: > 50 Ų)

### Abrindo no PyMOL

#### Opção 1: Usar o arquivo .pse (recomendado)

Se o PyMOL conseguiu gerar o arquivo .pse automaticamente:

```bash
pymol 1CRN_hydro.pse
```

A sessão abrirá com todas as representações e cores já configuradas!

#### Opção 2: Usar o script .pml

Se o arquivo .pse não foi gerado (PyMOL não estava no PATH durante a execução):

```bash
pymol 1CRN_hydro.pml
```

Ou abra o PyMOL e execute:

```python
# No console do PyMOL:
@1CRN_hydro.pml
```

#### Opção 3: Carregar manualmente o PDB

Abra o PyMOL e carregue o PDB anotado:

```python
# No console do PyMOL:
load 1CRN_hydro.pdb

# Aplicar gradiente de cores baseado no B-factor
spectrum b, blue_white_red, minimum=-4.5, maximum=4.5

# Mostrar representações
show cartoon
show surface
show sticks

# Ajustar transparência
set transparency, 0.5
```

### Representações Disponíveis

O script PyMOL gerado ativa automaticamente três representações:

1. **Cartoon (Fita)**: Mostra a estrutura secundária (α-hélices, β-folhas)
   - Colorida pelo gradiente do B-factor
   - Use `hide cartoon` para ocultar

2. **Sticks (Bastões)**: Mostra os átomos individuais
   - Colorido por tipo de átomo (C=cinza, N=azul, O=vermelho)
   - Use `hide sticks` para ocultar

3. **Surface (Superfície)**: Mostra a superfície molecular
   - Colorida pelo gradiente do B-factor
   - Transparência de 50% (ajustável)
   - Use `hide surface` para ocultar

### Comandos Úteis no PyMOL

**Ajustar visualização:**
```python
# Alterar transparência da superfície
set transparency, 0.7

# Ocultar/mostrar representações
hide surface
show surface

# Mudar esquema de cores
spectrum b, rainbow, minimum=-4.5, maximum=4.5
spectrum b, red_green_blue

# Fundo preto
bg_color black

# Renderizar imagem de alta qualidade
ray 2400, 2400
png imagem_hd.png, dpi=300
```

**Selecionar regiões específicas:**
```python
# Selecionar apenas resíduos hidrofóbicos (B-factor > 2.0)
select hydrophobic, b > 2.0
show spheres, hydrophobic

# Selecionar apenas resíduos hidrofílicos (B-factor < -2.0)
select hydrophilic, b < -2.0
show spheres, hydrophilic
```

### Comparando Hidrofobicidade e SASA

Você pode carregar ambas as análises simultaneamente:

```python
# Carregar hidrofobicidade
load 1CRN_hydro.pdb, hydro
spectrum b, blue_white_red, hydro, minimum=-4.5, maximum=4.5

# Carregar SASA
load 1CRN_sasa.pdb, sasa
spectrum b, blue_white_red, sasa, minimum=0, maximum=60

# Alinhar estruturas
align sasa, hydro

# Alternar entre visualizações
disable hydro
enable sasa
```

### Salvando Sessão Personalizada

Após ajustar a visualização no PyMOL, salve sua sessão:

```python
save minha_analise.pse
```

### Exportando Imagens

```python
# PNG de alta resolução
ray 2400, 2400
png figura1.png, dpi=300

# Renderização com ray tracing
set ray_shadow, 1
set ray_trace_mode, 1
ray

# Exportar como VRML (3D)
save estrutura.wrl
```

---

## Visualização com PyMOL

O BioHub pode gerar arquivos otimizados para visualização no PyMOL, incluindo PDBs anotados com propriedades no B-factor e scripts de visualização pré-configurados.

### Instalação do PyMOL

**macOS (via Homebrew):**
```bash
brew install pymol
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install pymol
```

**Windows:**
Baixe o instalador em [https://pymol.org/](https://pymol.org/)

**Verificar instalação:**
```bash
pymol -c
```

### Usando Arquivos Gerados pelo BioHub

#### Opção 1: Carregar PDB Anotado Manualmente

Após gerar um PDB com `--write-pdb`:

```bash
# Gerar PDB anotado com hidrofobicidade
python3 biohub.py exposure 1CRN.pdb --write-pdb 1CRN_hydro.pdb
```

Abrir no PyMOL:

```bash
pymol 1CRN_hydro.pdb
```

No console do PyMOL, aplicar o gradiente de cores:

```python
# Para hidrofobicidade (azul=hidrofílico, vermelho=hidrofóbico)
spectrum b, blue_white_red, minimum=-4.5, maximum=4.5

# Para SASA (azul=enterrado, vermelho=exposto)
spectrum b, blue_white_red, minimum=0, maximum=60

# Configurar visualização
show cartoon
show surface
set transparency, 0.3
```

#### Opção 2: Usar Script PyMOL (.pml)

O BioHub gera automaticamente um script `.pml` com visualização pré-configurada:

```bash
# Gerar PDB + script PyMOL
python3 biohub.py exposure 1CRN.pdb --write-pdb 1CRN_hydro.pdb --pymol 1CRN_hydro.pse
```

Arquivos gerados:
- `1CRN_hydro.pdb`: Estrutura com valores no B-factor
- `1CRN_hydro.pml`: Script com comandos de visualização
- `1CRN_hydro.pse`: Sessão PyMOL (se PyMOL estiver no PATH)

**Executar o script:**

```bash
pymol 1CRN_hydro.pml
```

O script configura automaticamente:
- Carregamento da estrutura
- Gradiente de cores apropriado
- Representação cartoon + surface
- Transparência e qualidade de renderização
- Centralização da visualização

#### Opção 3: Abrir Sessão PyMOL (.pse)

Se o PyMOL estiver instalado no PATH, o BioHub gera uma sessão completa:

```bash
pymol 1CRN_hydro.pse
```

A sessão abre com tudo pré-configurado, pronta para análise e geração de imagens.

### Exemplos Práticos

#### Visualizar Hidrofobicidade

```bash
# Gerar dados
python3 biohub.py exposure proteina.pdb --write-pdb proteina_hydro.pdb --pymol proteina_hydro.pse

# Abrir no PyMOL
pymol proteina_hydro.pml
```

**Personalizar no PyMOL:**
```python
# Ajustar transparência
set transparency, 0.5

# Mudar esquema de cores
spectrum b, rainbow, minimum=-4.5, maximum=4.5

# Destacar regiões hidrofóbicas (valores > 2.0)
select hydrophobic, b > 2.0
show spheres, hydrophobic
color red, hydrophobic

# Destacar regiões hidrofílicas (valores < -2.0)
select hydrophilic, b < -2.0
show spheres, hydrophilic
color blue, hydrophilic
```

#### Visualizar SASA

```bash
# Gerar dados (usar menos pontos para rapidez em demos)
python3 biohub.py sasa proteina.pdb --num-points 100 --write-pdb proteina_sasa.pdb --pymol proteina_sasa.pse

# Abrir no PyMOL
pymol proteina_sasa.pml
```

**Análise no PyMOL:**
```python
# Identificar resíduos completamente enterrados (SASA = 0)
select buried, b = 0
show sticks, buried
color gray, buried

# Identificar resíduos expostos (SASA > 30)
select exposed, b > 30
show surface, exposed
color red, exposed

# Gerar imagem de alta qualidade
set ray_trace_mode, 1
ray 2400, 2400
png proteina_sasa_rendered.png
```

#### Comparar Hidrofobicidade e SASA

```bash
# Gerar ambas as análises
python3 biohub.py exposure proteina.pdb --write-pdb proteina_hydro.pdb --pymol proteina_hydro.pse
python3 biohub.py sasa proteina.pdb --write-pdb proteina_sasa.pdb --pymol proteina_sasa.pse

# Abrir no PyMOL e carregar ambos
pymol proteina_hydro.pdb proteina_sasa.pdb
```

**No console do PyMOL:**
```python
# Configurar visualização lado a lado
set grid_mode, 1

# Objeto 1: Hidrofobicidade
spectrum b, blue_white_red, proteina_hydro, minimum=-4.5, maximum=4.5
disable proteina_hydro

# Objeto 2: SASA
spectrum b, blue_white_red, proteina_sasa, minimum=0, maximum=60
disable proteina_sasa

# Alternar entre visualizações
enable proteina_hydro  # Ver hidrofobicidade
disable proteina_hydro
enable proteina_sasa   # Ver SASA
```

### Salvar Sessão Manualmente

Se você fez modificações no PyMOL e quer salvar:

```python
# No console do PyMOL
save minha_sessao.pse
```

### Exportar Imagens

```python
# Renderização básica
png imagem.png, width=1200, height=1200

# Renderização de alta qualidade com ray tracing
set ray_trace_mode, 1
set antialias, 2
ray 2400, 2400
png imagem_hq.png
```

### Dicas de Visualização

**Esquemas de cores alternativos:**
```python
# Gradiente suave
spectrum b, rainbow

# Escala térmica
spectrum b, red_yellow_green_cyan_blue

# Preto e branco
spectrum b, black_white
```

**Representações úteis:**
```python
# Apenas superfície
hide everything
show surface

# Apenas cartoon
hide everything
show cartoon
set cartoon_fancy_helices, 1

# Combinação cartoon + sticks para resíduos chave
show cartoon
select key_residues, resi 10+25+45
show sticks, key_residues
```

**Análise quantitativa:**
```python
# Estatísticas do B-factor (valores calculados)
print "Min B-factor:", cmd.get_min("b")
print "Max B-factor:", cmd.get_max("b")
print "Mean B-factor:", cmd.get_mean("b")

# Contar átomos por faixa de valor
select low_hydro, b < -2
print "Átomos hidrofílicos:", cmd.count_atoms("low_hydro")

select high_hydro, b > 2
print "Átomos hidrofóbicos:", cmd.count_atoms("high_hydro")
```

---

## Estrutura de Saída dos Comandos

### fetchpdb
* **Terminal**: Mensagens de progresso e informações extraídas do cabeçalho PDB
* **Arquivo**: Arquivo PDB completo baixado do RCSB

### fasta
* **Terminal**: Sequência em formato FASTA (se `-o` não for especificado)
* **Arquivo**: Sequência em formato FASTA

### csv2fasta
* **Terminal**: Múltiplas sequências em formato FASTA (se `-o` não for especificado)
* **Arquivo**: Múltiplas sequências em formato FASTA

### physchem
* **Terminal**: Tabela formatada com propriedades e composição de aminoácidos
* **Arquivo CSV**: Duas colunas (Propriedade, Valor)

### contacts
* **Terminal**: Lista de contatos (Resíduo1, Resíduo2, Distância em Å)
* **Arquivo CSV**: Três colunas (Residuo1, Residuo2, Distancia_A)

### exposure
* **Terminal**: Tabela com Chain, ResNum, ResName, AtomNum, AtomName e valores de hidrofobicidade
* **Arquivo CSV**: Seis colunas (Chain, ResNum, ResName, AtomNum, AtomName, Hydrophobicity)
* **Arquivo PDB** (com `--write-pdb`): Estrutura 3D com hidrofobicidade no B-factor

### sasa
* **Terminal**: SASA total da molécula + tabela com SASA por átomo
* **Arquivo CSV**: Seis colunas (Chain, ResNum, ResName, AtomNum, AtomName, SASA_A2)
* **Arquivo PDB** (com `--write-pdb`): Estrutura 3D com SASA no B-factor

### apbs [BETA, ainda em validação]
* **Terminal**: Progresso da execução e energia de solvatação final em kJ/mol
* **Arquivos** (com `--no-cleanup`): PQR, apbs.in, mapas de potencial eletrostático

---

## Expansão Futura

Está em planejamento o desenvolvimento de uma versão web da BioHub. A aplicação web terá uma interface gráfica intuitiva para realizar as mesmas análises, com visualizações interativas dos resultados, tornando os outputs acessíveis a um público mais amplo. A lógica de cálculo será portada para JavaScript, permitindo análises rápidas diretamente no navegador.

Uma versão beta está disponível em: [(https://madsondeluna.github.io/biohub-beta/)](https://madsondeluna.github.io/apps/biohub/index.html)

---

## Licença

Este projeto está licenciado sob a **Licença MIT**. Veja o arquivo `LICENSE` para mais detalhes.

---

```
      |\      _,,,---,,_
ZZZzz /,`.-'`'    -.  ;-;;,_
     |,4-  ) )-,_. ,\ (  `'-'
    '---''(_/--'  `-'\_)
```

---

Made with ❤️, a shared caffeine addiction and collective chaos from the BioHub WhatsApp group.
