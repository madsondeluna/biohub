# BioHub: Uma Plataforma para Análise de Sequências e Estruturas de Proteínas

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.x-brightgreen.svg)](https://www.python.org/)

**BioHub** é uma ferramenta de linha de comando leve e livre de dependências, escrita em Python, para realizar análises bioinformáticas a partir de arquivos de estrutura de proteínas (PDB) e sequências de aminoácidos.

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

## Equipe Desenvolvedora (Bioinformática - UFMG)

- Alisson Clementino da Silva 
- Ana Karoline da Nóbrega Nunes Alves
- Laura da Silva Ribeiro de Souza
- Leonardo Henrique da Silva
- Madson Allan de Luna Aragão

---

## Visão Geral

Este projeto foi desenvolvido com o objetivo de fornecer uma solução simples, rápida e altamente portátil para tarefas comuns de análise de proteínas. É ideal para:

* **Fins educacionais**: Demonstra a implementação de algoritmos bioinformáticos fundamentais sem a abstração de bibliotecas complexas.
* **Integração em pipelines**: A natureza leve e a interface de linha de comando facilitam a automação de análises em scripts de shell.
* **Acessibilidade**: Por não requerer bibliotecas externas, a ferramenta pode ser executada em qualquer ambiente com uma instalação padrão do Python 3, eliminando problemas de gerenciamento de dependências.

---

## Como Funciona

A BioHub foi intencionalmente construída utilizando apenas a biblioteca padrão do Python. Constantes físico-químicas, como pesos moleculares e escalas de hidropatia, são armazenadas localmente em dicionários.

**Justificativa**: Esta abordagem garante a máxima portabilidade e simplicidade. Embora o uso de bancos de dados externos ou pacotes especializados pudesse oferecer maior precisão ou dados mais atualizados, a implementação atual é suficiente para fins demonstrativos e educacionais, focando na lógica algorítmica.

---

## Funcionalidades Detalhadas

A ferramenta é organizada em quatro subcomandos principais:

### 1. `fasta`

Converte um arquivo PDB em uma sequência no formato FASTA. O processo foca nos registros `ATOM` e extrai a sequência de aminoácidos da primeira cadeia (chain) encontrada no arquivo.

### 2. `physchem`

Calcula um conjunto de propriedades físico-químicas essenciais para uma dada sequência de aminoácidos:

* **Peso Molecular (MW)**: A massa da molécula, fundamental para experimentos de espectrometria de massa e SDS-PAGE.
* **Ponto Isoelétrico (pI)**: O pH no qual a proteína tem carga líquida zero, influenciando sua solubilidade e comportamento em cromatografia de troca iônica.
* **Coeficiente de Extinção**: Estima como a proteína absorve luz a 280 nm, útil para determinar a concentração da proteína via espectrofotometria.
* **GRAVY (Grand Average of Hydropathicity)**: Um índice do caráter hidrofóbico ou hidrofílico geral de uma proteína. Valores positivos indicam hidrofobicidade (ex.: proteínas de membrana), enquanto valores negativos indicam hidrofilicidade (ex.: proteínas citosólicas).
* **Composição de Aminoácidos**: A contagem e a frequência de cada resíduo.

### 3. `contacts`

Identifica e lista contatos intramoleculares com base na distância entre os átomos de **Carbono Alfa (CA)**. A distância entre os CAs é um excelente proxy para a proximidade entre os resíduos. Este cálculo é útil para:

* Identificar o núcleo hidrofóbico (core) da proteína.
* Analisar a topologia do enovelamento.
* Estudar interações de longo alcance que estabilizam a estrutura terciária.

### 4. `exposure`

Prevê regiões potencialmente expostas ao solvente ou enterradas no interior da proteína. Utiliza a **escala de hidropatia de Kyte-Doolittle** com um método de janela deslizante. Para cada resíduo, calcula-se a média de hidropatia dos resíduos vizinhos. Pontuações altas indicam regiões hidrofóbicas (provavelmente internas), enquanto pontuações baixas indicam regiões hidrofílicas (provavelmente na superfície).

### 5. `sasa`
Calcula a Área de Superfície Acessível ao Solvente (SASA). Esta funcionalidade implementa o algoritmo de Shrake-Rupley em Python puro para estimar a área da superfície da proteína que está em contato com o solvente. É uma métrica fundamental para estudos de enovelamento, estabilidade e interações moleculares.

### 6. `apbs`
Calcula a energia de solvatação eletrostática. Este comando atua como um wrapper, automatizando a execução dos softwares PDB2PQR e APBS. Ele prepara os arquivos necessários, executa os cálculos de Poisson-Boltzmann e extrai o valor final da energia, simplificando uma análise computacionalmente complexa.

---

## Requisitos

* Python 3.x

Nenhuma outra biblioteca é necessária.

---

## Instalação

1. Clone este repositório:

   ```bash
   git clone https://github.com/seu-usuario/biohub.git
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
python biohub.py [COMANDO] [ARGUMENTOS] [OPÇÕES]
```

Ou, se o arquivo for executável:

```bash
./biohub.py [COMANDO] [ARGUMENTOS] [OPÇÕES]
```

### Ajuda (ou `biohub.py -h`)

```bash
usage: biohub.py [-h] {fasta,physchem,contacts,exposure} ...

Ferramenta de bioinformática para análise de PDB e sequências.

positional arguments:
  {fasta,physchem,contacts,exposure}
                        Função a ser executada
    fasta               Converte um arquivo PDB em uma sequência FASTA.
    physchem            Calcula propriedades físico-químicas de uma
                        sequência.
    contacts            Calcula contatos intramoleculares a partir de um
                        arquivo PDB.
    exposure            Prevê a exposição ao solvente com a escala Kyte-
                        Doolittle.

optional arguments:
  -h, --help            show this help message and exit
```

### Exemplos Detalhados

Suponha que você tenha um arquivo `proteina.pdb`.

#### 1) Converter PDB para FASTA

Exibir a sequência no terminal:

```bash
python biohub.py fasta proteina.pdb
```

Salvar a sequência em um arquivo (`-o` ou `--output`):

```bash
python biohub.py fasta proteina.pdb -o proteina.fasta
```

#### 2) Calcular Propriedades Físico-Químicas

Analisar uma sequência fornecida diretamente:

```bash
python biohub.py physchem "MKTAYIAKQ"
```

#### 3) Calcular Contatos Intramoleculares

Encontrar contatos com o limiar padrão (8.0 Å):

```bash
python biohub.py contacts proteina.pdb
```

Ajustar o limiar de distância para 10.5 Å (`-t` ou `--threshold`):

```bash
python biohub.py contacts proteina.pdb -t 10.5
```

#### 4) Prever Exposição ao Solvente

Executar com a janela padrão (9 resíduos):

```bash
python biohub.py exposure proteina.pdb
```

Usar uma janela maior (ex.: 19) para suavizar o perfil e identificar regiões extensas (ex.: hélices transmembrana) com `-w` ou `--window`:

```bash
python biohub.py exposure proteina.pdb -w 19
```

#### 5) Calcular SASA (Área de Superfície Acessível ao Solvente)

Executar com parâmetros padrão (sonda de 1.4 Å, 960 pontos):

```bash
python biohub.py sasa proteina.pdb
```

Ajustar a precisão e o raio da sonda:

```bash
python biohub.py sasa proteina.pdb --probe-radius 1.5 --num-points 2000
```

#### 6) Calcular Energia de Solvatação (APBS)

Executar a análise padrão (arquivos temporários são removidos):

```bash
python biohub.py apbs proteina.pdb
```

Executar e manter os arquivos intermediários (PQR, apbs.in, etc.) para depuração:

```bash
python biohub.py apbs proteina.pdb --no-cleanup
```


---

## Automação com Scripts

A BioHub é ideal para uso em scripts de shell para processar múltiplos arquivos em lote.

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
    python biohub.py contacts "$pdb" > "${OUTPUT_DIR}/${base_name}_contacts.txt"

    # Gera o perfil de exposição ao solvente
    python biohub.py exposure "$pdb" > "${OUTPUT_DIR}/${base_name}_exposure.txt"

done

echo "Análise concluída. Resultados em ${OUTPUT_DIR}/"
```

---

## Expansão Futura

Está em planejamento o desenvolvimento de uma versão web da BioHub. A aplicação web terá uma interface gráfica intuitiva para realizar as mesmas análises, com visualizações interativas dos resultados, tornando-os acessíveis a um público mais amplo. A lógica de cálculo será portada para JavaScript, permitindo análises rápidas diretamente no navegador.

Uma versão beta está disponível em: https://madsondeluna.github.io/biohub-beta/

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
