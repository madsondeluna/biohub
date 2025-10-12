# BioHub: Uma Plataforma para Análise de Sequências e Estruturas de Proteínas

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.x-brightgreen.svg)](https://www.python.org/)
[![CLI](https://img.shields.io/badge/Interface-CLI-orange.svg)](#)
[![PDB](https://img.shields.io/badge/Format-PDB-blue)](https://www.wwpdb.org/)
[![CSV](https://img.shields.io/badge/Format-CSV-yellow)](https://en.wikipedia.org/wiki/Comma-separated_values)
[![FASTA](https://img.shields.io/badge/Format-FASTA-lightgrey)](https://en.wikipedia.org/wiki/FASTA_format)


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

Converte arquivos CSV contendo sequências de aminoácidos para o formato FASTA. Esta funcionalidade é especialmente útil para processar grandes volumes de sequências armazenadas em planilhas ou bancos de dados tabulares.

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

Prevê regiões potencialmente expostas ao solvente ou enterradas no interior da proteína. Utiliza a **escala de hidropatia de Kyte-Doolittle** com um método de janela deslizante. Para cada resíduo, calcula-se a média de hidropatia dos resíduos vizinhos. Pontuações altas indicam regiões hidrofóbicas (provavelmente internas), enquanto pontuações baixas indicam regiões hidrofílicas (provavelmente na superfície).

### 5. `sasa`
Calcula a Área de Superfície Acessível ao Solvente (SASA). Esta funcionalidade implementa o algoritmo de Shrake-Rupley para estimar a área da superfície da proteína que está em contato com o solvente. É uma métrica fundamental para estudos de enovelamento, estabilidade e interações moleculares.

### 7. `sasa`

Calcula a **Área de Superfície Acessível ao Solvente (SASA)** usando o algoritmo de Shrake-Rupley implementado em Python puro. Esta métrica é fundamental para estudos de enovelamento, estabilidade e interações moleculares.

**Características do cálculo:**
* Raio da sonda customizável (padrão: 1.4 Å para água)
* Número de pontos ajustável para controle de precisão (padrão: 960 pontos)
* Cálculo por resíduo e SASA total da molécula
* Uso de raios de Van der Waals específicos para cada tipo de átomo

### 8. `apbs`

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
usage: biohub.py [-h] {fetchpdb,fasta,csv2fasta,physchem,contacts,exposure,sasa,apbs} ...

BioHub: Uma ferramenta CLI para análise de proteínas.

positional arguments:
  {fetchpdb,fasta,csv2fasta,physchem,contacts,exposure,sasa,apbs}
                        Função a ser executada
    fetchpdb            Baixa um arquivo PDB do RCSB.
    fasta               Converte um arquivo PDB em uma sequência FASTA.
    csv2fasta           Converte um arquivo CSV em um formato FASTA.
    physchem            Calcula propriedades físico-químicas de uma sequência.
    contacts            Calcula contatos intramoleculares a partir de um arquivo PDB.
    exposure            Prevê a exposição ao solvente com a escala Kyte-Doolittle.
    sasa                Calcula a Área de Superfície Acessível ao Solvente (SASA).
    apbs                Calcula a energia de solvatação eletrostática.

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

**Ajuda do comando:**
```bash
python3 biohub.py fetchpdb -h
```

**Opções disponíveis:**
* `PDBID`: O código de 4 caracteres do PDB (obrigatório)
* `-o, --output ARQUIVO`: Nome do arquivo de saída (padrão: `PDBID.pdb`)

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

### 6) Prever Exposição ao Solvente

Executar com a janela padrão (9 resíduos):

```bash
python3 biohub.py exposure proteina.pdb
```

Usar uma janela maior (19 resíduos) para suavizar o perfil:

```bash
python3 biohub.py exposure proteina.pdb -w 19
```

Salvar resultados em CSV:

```bash
python3 biohub.py exposure proteina.pdb -w 9 -o exposicao.csv
```

**Ajuda do comando:**
```bash
python3 biohub.py exposure -h
```

**Opções disponíveis:**
* `ARQUIVO_PDB`: Caminho para o arquivo PDB de entrada (obrigatório)
* `-w, --window INT`: Tamanho da janela deslizante (deve ser ímpar, padrão: 9)
* `-o, --output ARQUIVO_CSV`: Salva os resultados em um arquivo CSV

**Interpretação dos resultados:**
* **Scores positivos**: Regiões hidrofóbicas (provavelmente enterradas)
* **Scores negativos**: Regiões hidrofílicas (provavelmente expostas ao solvente)

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

Salvar resultados por resíduo em CSV:

```bash
python3 biohub.py sasa proteina.pdb -o sasa_resultados.csv
```

**Ajuda do comando:**
```bash
python3 biohub.py sasa -h
```

**Opções disponíveis:**
* `ARQUIVO_PDB`: Caminho para o arquivo PDB de entrada (obrigatório)
* `--probe-radius FLOAT`: Raio da sonda do solvente em Angstroms (padrão: 1.4 para água)
* `--num-points INT`: Número de pontos na superfície de cada átomo (padrão: 960)
* `-o, --output ARQUIVO_CSV`: Salva os resultados por resíduo em um arquivo CSV

**Observações:**
* Valores mais altos de `--num-points` aumentam a precisão mas também o tempo de processamento
* O raio de 1.4 Å é o padrão para moléculas de água
* A SASA total da molécula é sempre exibida no stderr

---

### 8) Calcular Energia de Solvatação (APBS)

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

## Automação com Scripts

A BioHub é ideal para uso em scripts de shell para processar múltiplos arquivos em lote.

### Exemplo 1: Análise em Lote de Múltiplos PDBs

```bash
#!/bin/bash

# Define o diretório de saída
OUTPUT_DIR="analise_resultados"
mkdir -p "$OUTPUT_DIR"

# Itera sobre todos os arquivos PDB no diretório atual
for pdb in *.pdb; do
    base_name=$(basename "$pdb" .pdb)
    echo "Processando ${base_name}..."

    # Gera o arquivo de contatos
    python3 biohub.py contacts "$pdb" -o "${OUTPUT_DIR}/${base_name}_contacts.csv"

    # Gera o perfil de exposição ao solvente
    python3 biohub.py exposure "$pdb" -o "${OUTPUT_DIR}/${base_name}_exposure.csv"

    # Calcula SASA
    python3 biohub.py sasa "$pdb" -o "${OUTPUT_DIR}/${base_name}_sasa.csv"

done

echo "Análise concluída. Resultados em ${OUTPUT_DIR}/"
```

### Exemplo 2: Download e Análise Automatizada

```bash
#!/bin/bash

# Lista de códigos PDB para analisar
PDB_IDS=("1A2B" "2GHI" "3JKL" "4MNO")

for pdb_id in "${PDB_IDS[@]}"; do
    echo "Processando ${pdb_id}..."

    # Baixa o arquivo PDB
    python3 biohub.py fetchpdb "$pdb_id"

    # Converte para FASTA
    python3 biohub.py fasta "${pdb_id}.pdb" -o "${pdb_id}.fasta"

    # Calcula propriedades físico-químicas
    sequence=$(tail -n 1 "${pdb_id}.fasta")
    python3 biohub.py physchem "$sequence" -o "${pdb_id}_physchem.csv"

    echo "Concluído: ${pdb_id}"
done
```

### Exemplo 3: Processamento de Sequências de CSV

```bash
#!/bin/bash

# Converte CSV para FASTA
python3 biohub.py csv2fasta sequencias.csv --header --id-col "ProteinID" --seq-col "Sequence" -o todas_sequencias.fasta

# Extrai cada sequência e calcula propriedades
grep -A1 "^>" todas_sequencias.fasta | while read header; do
    read sequence
    if [[ $header == ">"* ]]; then
        id=$(echo $header | sed 's/>//')
        python3 biohub.py physchem "$sequence" -o "${id}_properties.csv"
    fi
done
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
* **Terminal**: Tabela com posição, resíduo e score de hidropatia
* **Arquivo CSV**: Três colunas (Posicao, Residuo, Score_Hidropatia)

### sasa
* **Terminal**: SASA total + tabela com SASA por resíduo
* **Arquivo CSV**: Duas colunas (Residuo, SASA_A2)

### apbs
* **Terminal**: Progresso da execução e energia de solvatação final em kJ/mol
* **Arquivos** (com `--no-cleanup`): PQR, apbs.in, mapas de potencial eletrostático

---

## Expansão Futura

Está em planejamento o desenvolvimento de uma versão web da BioHub. A aplicação web terá uma interface gráfica intuitiva para realizar as mesmas análises, com visualizações interativas dos resultados, tornando-os acessíveis a um público mais amplo. A lógica de cálculo será portada para JavaScript, permitindo análises rápidas diretamente no navegador.

Uma versão beta está disponível em: [https://madsondeluna.github.io/biohub-beta/](https://madsondeluna.github.io/biohub-beta/)

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

***Made with ❤️ and countless cups of coffee!***
