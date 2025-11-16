# BioHub - Guia Detalhado de Funções para Apresentação

Este documento organiza todas as funções do BioHub ([biohub.py](../biohub.py) e [biohub_viz.py](../biohub_viz.py)) de forma didática para apresentação, separando:

1. **Dicionários e Constantes** - Dados científicos necessários
2. **Função de Cálculo** - Implementação do algoritmo
3. **Formatação de Output** - Como os resultados são apresentados
4. **Argumentos CLI** - Parâmetros que o usuário pode usar
5. **Detalhes de Implementação** - Como cada flag/argumento funciona internamente

## Resumo Executivo

### Visao Geral das Funcoes

| Funcao | Proposito | Flags Principais | Codigo Chamado |
|--------|-----------|------------------|----------------|
| fetchpdb | Download de estruturas PDB | --chains, --protein-only | handle_fetch_pdb() + filter_pdb_content() |
| fasta | Conversao PDB para FASTA | -o | get_sequence_from_pdb() |
| csv2fasta | Conversao CSV para FASTA | --header, --delimiter | handle_csv_to_fasta() |
| physchem | Propriedades fisico-quimicas | --plot-treemap, --plot-hydro | calculate_physicochemical_properties() |
| contacts | Analise de contatos | -t, --plot | calculate_intramolecular_contacts() |
| hydrophoby | Perfil de hidrofobicidade | --write-pdb, --pymol | predict_solvent_hydrophoby() |
| sasa | Superficie acessivel ao solvente | --probe-radius, --num-points | calculate_sasa() (Shrake-Rupley) |

### Fluxo Tipico de Uso

```bash
# 1. Download da estrutura PDB, filtrando cadeia A e removendo agua
python biohub.py fetchpdb 1TUP --chains A --protein-only -o 1TUP_clean.pdb

# 2. Extrair sequencia em formato FASTA
python biohub.py fasta 1TUP_clean.pdb -o 1TUP.fasta

# 3. Calcular propriedades fisico-quimicas da sequencia
python biohub.py physchem $(grep -v ">" 1TUP.fasta) -o propriedades.csv

# 4. Calcular SASA e gerar visualizacao PyMOL
python biohub.py sasa 1TUP_clean.pdb --write-pdb 1TUP_sasa.pdb --pymol sasa.pse

# 5. Calcular hidrofobicidade e gerar perfil
python biohub.py hydrophoby 1TUP_clean.pdb --write-pdb 1TUP_hydro.pdb --plot-hydrophoby perfil.png

# 6. Analisar contatos intramoleculares
python biohub.py contacts 1TUP_clean.pdb -t 8.0 --plot contact_map.png
```

---

## 00 - fetchpdb: Download Automático de Estruturas PDB

### 1. Dicionários e Constantes

Nenhum dicionário específico. Usa constantes globais:
- **URL base RCSB**: `https://files.rcsb.org/download/{PDB_ID}.pdb`

### 2. Função de Cálculo

**Arquivo**: [biohub.py:269-316](../biohub.py#L269-L316)

```python
def handle_fetch_pdb(args):
    """Baixa um arquivo PDB do banco de dados RCSB PDB."""
    pdb_id = args.pdb_id.upper()

    # Validação do ID (4 caracteres)
    if len(pdb_id) != 4:
        print(f"Erro: O ID do PDB '{pdb_id}' é inválido.")
        return

    # Monta URL e faz download
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_file = args.output if args.output else f"{pdb_id}.pdb"

    with urllib.request.urlopen(url) as response:
        data = response.read()
        # Salva arquivo

    # Aplica filtros opcionais (cadeias, protein-only)
    if chains_list or args.protein_only:
        atoms_kept, atoms_removed = filter_pdb_content(...)
```

**Função auxiliar - Filtros**: [biohub.py:212-267](../biohub.py#L212-L267)
```python
def filter_pdb_content(pdb_filepath, chains=None, protein_only=False):
    """Filtra PDB removendo água, ligantes ou selecionando cadeias."""
    # Remove HETATM se protein_only
    # Filtra por chain_id se especificado
    # Mantém apenas linhas ATOM desejadas
```

### 3. Formatação de Output

**Output primário**: Arquivo PDB salvo
**Output secundário**: Informações do cabeçalho extraídas ([biohub.py:162-184](../biohub.py#L162-L184))

```python
def extract_pdb_header_info(pdb_filepath):
    """Extrai: classification, dep_date, title, organism."""
    # Lê linhas HEADER, TITLE, SOURCE do PDB
    # Retorna dicionário com metadados
```

Exibição no terminal:
```
--- Informações da Estrutura ---
  Título    : TUMOR SUPPRESSOR PROTEIN P53
  Data      : 24-JAN-95
  Classe    : TUMOR SUPPRESSOR
  Organismo : HOMO SAPIENS
```

### 4. Argumentos CLI

**Comando**: `biohub.py fetchpdb`

| Argumento | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `pdb_id` | Posicional | - | ID de 4 caracteres (ex: 1TUP) |
| `-o, --output` | Opcional | `{PDB_ID}.pdb` | Nome do arquivo de saída |
| `--chains` | Opcional | Todas | Cadeias a manter (ex: A,B) |
| `--protein-only` | Flag | False | Remove água e ligantes |

**Exemplo de uso**:
```bash
python biohub.py fetchpdb 1TUP -o proteina.pdb --chains A --protein-only
```

### 5. Detalhes de Implementacao de Cada Flag

#### Flag: `--chains`
**Arquivo**: [biohub.py:1006](../biohub.py#L1006)
**Parser**: [biohub.py:289-290](../biohub.py#L289-L290)
**Funcao chamada**: `filter_pdb_content()` em [biohub.py:212-267](../biohub.py#L212-L267)

**Como funciona**:
1. O argumento eh parseado na linha 289: `chains_list = [c.strip().upper() for c in args.chains.split(',')]`
   - Divide a string por virgula (ex: "A,B" vira ['A', 'B'])
   - Remove espacos em branco (.strip())
   - Converte para maiusculo (.upper())

2. A funcao `filter_pdb_content()` recebe a lista de cadeias e:
   - Le o arquivo PDB linha por linha (linha 228)
   - Para cada linha ATOM ou HETATM (linha 233):
     - Extrai o chain_id da coluna 21 do formato PDB (linha 246)
     - Verifica se o chain_id esta na lista de cadeias desejadas (linha 247)
     - Se NAO estiver, marca keep_line = False (linha 248)
     - Incrementa contador atoms_removed (linha 249)

3. Apenas linhas com keep_line = True sao mantidas (linha 260-261)
4. O arquivo eh reescrito com o conteudo filtrado (linha 264-265)

**Formato PDB relevante**:
```
ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N
                    ^
                    Coluna 21 = Chain ID
```

#### Flag: `--protein-only`
**Arquivo**: [biohub.py:1007](../biohub.py#L1007)
**Funcao chamada**: `filter_pdb_content()` em [biohub.py:212-267](../biohub.py#L212-L267)

**Como funciona**:
1. Quando protein_only=True, a funcao aplica dois filtros:

   **Filtro 1 - Remove HETATM** (linhas 235-237):
   - Se a linha comeca com "HETATM", marca keep_line = False
   - HETATM = heteroatomos (ligantes, ions, moleculas pequenas)
   - ATOM = atomos da cadeia polipeptidica principal

   **Filtro 2 - Remove agua (HOH)** (linhas 240-242):
   - Extrai o nome do residuo das colunas 17-20 do formato PDB
   - Se for "HOH" (agua), marca keep_line = False
   - HOH = codigo PDB para molecula de agua

2. Tambem remove linhas HETATM do cabecalho (linhas 255-258)

**Exemplo de linhas removidas**:
```
HETATM 3046  O   HOH A 301      15.234  28.123  10.456  1.00 30.12           O  <- Agua (removida)
HETATM 3047 MG    MG A 302      25.678  32.901   8.234  1.00 25.67          MG  <- Ion magnesio (removido)
HETATM 3048  C1  ATP A 303      18.456  35.234  12.567  1.00 40.23           C  <- Ligante ATP (removido)
ATOM   3049  N   ALA A  45      22.345  30.123   9.876  1.00 35.45           N  <- Proteina (MANTIDA)
```

#### Flag: `-o, --output`
**Arquivo**: [biohub.py:1005](../biohub.py#L1005)
**Uso**: [biohub.py:277](../biohub.py#L277)

**Como funciona**:
- Linha 277: `output_file = args.output if args.output else f"{pdb_id}.pdb"`
- Se o usuario NAO especificar -o, usa o padrao: {PDB_ID}.pdb
- Se o usuario especificar, usa o nome fornecido
- Exemplo: `fetchpdb 1TUP` salva como "1TUP.pdb"
- Exemplo: `fetchpdb 1TUP -o minha_proteina.pdb` salva como "minha_proteina.pdb"

#### Argumento: `pdb_id`
**Arquivo**: [biohub.py:1004](../biohub.py#L1004)
**Validacao**: [biohub.py:271-275](../biohub.py#L271-L275)
**Download**: [biohub.py:276-284](../biohub.py#L276-L284)

**Como funciona**:
1. Converte para maiusculo (linha 271): `pdb_id = args.pdb_id.upper()`
2. Valida tamanho (linhas 273-275):
   - Deve ter exatamente 4 caracteres
   - Se invalido, exibe erro e retorna

3. Monta URL do RCSB PDB (linha 276):
   - Template: `https://files.rcsb.org/download/{PDB_ID}.pdb`
   - Exemplo: `https://files.rcsb.org/download/1TUP.pdb`

4. Faz download usando urllib (linhas 281-284):
   - Abre conexao HTTP com urllib.request.urlopen()
   - Le todos os bytes (response.read())
   - Verifica se arquivo nao esta vazio
   - Escreve no arquivo de saida em modo binario ('wb')

5. Trata erros (linhas 311-316):
   - HTTPError 404: PDB ID nao encontrado
   - Remove arquivo vazio se download falhar

---

## 01 - fasta: Conversão PDB → FASTA

### 1. Dicionários e Constantes

**THREE_TO_ONE**: [biohub.py:31-36](../biohub.py#L31-L36)

```python
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
```

### 2. Função de Cálculo

**Arquivo**: [biohub.py:136-160](../biohub.py#L136-L160)

```python
def get_sequence_from_pdb(pdb_filepath: str) -> str:
    """Extrai a sequência de aminoácidos da primeira cadeia."""
    sequence = ""
    processed_residues = set()  # Evita duplicatas
    first_chain_id = None

    with open(pdb_filepath, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                # Define primeira cadeia encontrada
                if first_chain_id is None:
                    first_chain_id = line[21]

                # Processa apenas primeira cadeia
                if line[21] == first_chain_id:
                    res_id = (int(line[22:26]), line[21])
                    if res_id not in processed_residues:
                        res_name = line[17:20]
                        # Converte 3 letras → 1 letra
                        if res_name in THREE_TO_ONE:
                            sequence += THREE_TO_ONE[res_name]
                            processed_residues.add(res_id)
    return sequence
```

### 3. Formatação de Output

**Handler**: [biohub.py:318-330](../biohub.py#L318-L330)

```python
def handle_pdb_to_fasta(args):
    """Converte PDB para formato FASTA."""
    sequence = get_sequence_from_pdb(args.pdb_file)
    header = f">sequence_from_{os.path.basename(args.pdb_file)}"
    fasta_output = f"{header}\n{sequence}\n"

    # Salva em arquivo ou stdout
    if args.output:
        with open(args.output, 'w') as f:
            f.write(fasta_output)
```

**Formato FASTA**:
```
>sequence_from_1TUP.pdb
MTAMEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDLLLF
```

### 4. Argumentos CLI

**Comando**: `biohub.py fasta`

| Argumento | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `pdb_file` | Posicional | - | Arquivo PDB de entrada |
| `-o, --output` | Opcional | stdout | Arquivo FASTA de saída |

**Exemplo de uso**:
```bash
python biohub.py fasta 1TUP.pdb -o sequencia.fasta
```

### 5. Detalhes de Implementacao

#### Argumento: `pdb_file`

**Parser CLI**: [biohub.py:1011](../biohub.py#L1011)
**Handler**: [biohub.py:318-330](../biohub.py#L318-L330)
**Funcao de extracao**: `get_sequence_from_pdb()` em [biohub.py:136-160](../biohub.py#L136-L160)

**Como funciona**:

1. A funcao `get_sequence_from_pdb()` le o arquivo PDB linha por linha (linha 143)

2. Para cada linha que comeca com "ATOM" (linha 144):
   - Extrai chain_id da coluna 21 (linha 145)
   - Na primeira linha ATOM, define first_chain_id (linha 147)
   - **Apenas processa atomos da primeira cadeia encontrada** (linha 150)

3. Evita duplicacao de residuos (linhas 139, 152):
   - Cria um set() chamado processed_residues
   - Usa tupla (res_num, chain_id) como identificador unico
   - Cada residuo tem varios atomos (N, CA, C, O, CB, etc.)
   - Adiciona o residuo apenas uma vez a sequencia

4. Conversao 3 letras para 1 letra (linhas 153-157):
   - Extrai res_name das colunas 17-20 (ex: "ALA")
   - Busca no dicionario THREE_TO_ONE (ex: "ALA" -> "A")
   - Concatena na string sequence

5. Formato PDB relevante:
```
ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N
ATOM      2  CA  MET A   1      21.489  30.263   5.552  1.00 48.12           C
ATOM      3  C   MET A   1      22.374  29.134   6.089  1.00 47.89           C
             ^      ^   ^   ^
          Col 17  Col 21 Col 22-26
          (res)  (chain) (res_num)
```
   - Colunas 17-20: Nome do residuo (MET)
   - Coluna 21: Chain ID (A)
   - Colunas 22-26: Numero do residuo (1)

**Por que apenas a primeira cadeia?**
- Arquivos PDB podem ter multiplas cadeias (A, B, C, etc.)
- Evita sequencias duplicadas em estruturas oligomericas
- Usuario pode filtrar cadeias com `fetchpdb --chains A` antes

#### Flag: `-o, --output`

**Parser CLI**: [biohub.py:1012](../biohub.py#L1012)
**Uso**: [biohub.py:326-330](../biohub.py#L326-L330)

**Como funciona**:

1. Se args.output existe (linha 326):
   - Abre arquivo em modo escrita (linha 327)
   - Escreve cabecalho FASTA + sequencia
   - Exibe mensagem de confirmacao no stderr (linha 328)

2. Se args.output nao existe (linha 329-330):
   - Imprime diretamente no stdout
   - Permite uso com pipes: `python biohub.py fasta 1TUP.pdb | grep ">"``

3. Formato do cabecalho FASTA (linha 323):
   - Template: `>sequence_from_{nome_do_arquivo}`
   - Usa os.path.basename() para extrair apenas o nome do arquivo
   - Exemplo: `>sequence_from_1TUP.pdb`

---

## 02 - physchem: Propriedades Físico-Químicas Completas

### 1. Dicionários e Constantes

#### MOLECULAR_WEIGHT: [biohub.py:38-43](../biohub.py#L38-L43)
```python
MOLECULAR_WEIGHT = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10,
    'C': 121.16, 'E': 147.13, 'Q': 146.15, 'G': 75.07,
    # ... (em Daltons)
}
```

#### KYTE_DOOLITTLE: [biohub.py:44-49](../biohub.py#L44-L49)
```python
KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5,
    'C': 2.5, 'E': -3.5, 'Q': -3.5, 'G': -0.4,
    # ... (valores positivos = hidrofóbico)
}
```

#### PKA_VALUES: [biohub.py:50-54](../biohub.py#L50-L54)
```python
PKA_VALUES = {
    'C-term': 3.65, 'N-term': 8.0,
    'D': 3.9, 'E': 4.07, 'H': 6.5,
    'C': 8.5, 'Y': 10.0, 'K': 10.0, 'R': 12.0
}
```

#### DIWV (Matriz de Instabilidade): [biohub.py:60-82](../biohub.py#L60-L82)
```python
DIWV = {
    'A': {'A': 1.0, 'C': 44.94, 'E': 1.0, ...},
    'C': {'A': 1.0, 'C': 1.0, 'E': 1.0, ...},
    # ... matriz 20x20 de pesos dipeptídicos
}
```

### 2. Função de Cálculo

**Arquivo**: [biohub.py:371-436](../biohub.py#L371-L436)

```python
def calculate_physicochemical_properties(args):
    """Calcula propriedades físico-químicas da sequência."""
    sequence = args.sequence.upper()
    length = len(sequence)

    # 1. Composição de aminoácidos
    aa_composition = {aa: sequence.count(aa) for aa in MOLECULAR_WEIGHT.keys()}

    # 2. Peso Molecular (Da)
    # Soma dos pesos - água perdida nas ligações peptídicas
    mw = sum(MOLECULAR_WEIGHT.get(aa, 0) * count
             for aa, count in aa_composition.items()) - (length - 1) * 18.015

    # 3. GRAVY (Hidropaticidade média)
    gravy = sum(KYTE_DOOLITTLE.get(aa, 0) for aa in sequence) / length

    # 4. Ponto Isoelétrico (pI)
    # Busca binária por pH onde carga líquida ≈ 0
    for ph in range(0, 1401):  # pH 0.00 a 14.00
        ph = ph * 0.01
        # Carga do N-terminal (positiva)
        net_charge = (10**PKA_VALUES['N-term']) / (10**PKA_VALUES['N-term'] + 10**ph)
        # Carga do C-terminal (negativa)
        net_charge -= (10**ph) / (10**PKA_VALUES['C-term'] + 10**ph)
        # Cadeias laterais ionizáveis
        for aa, count in aa_composition.items():
            if aa in ['R', 'H', 'K']:  # Básicos (positivos)
                net_charge += count * (10**PKA_VALUES[aa]) / (10**PKA_VALUES[aa] + 10**ph)
            elif aa in ['D', 'E', 'C', 'Y']:  # Ácidos (negativos)
                net_charge -= count * (10**ph) / (10**PKA_VALUES[aa] + 10**ph)
        # Guarda pH com carga mais próxima de zero
        if abs(net_charge) < min_charge_diff:
            pi = ph

    # 5. Índice Alifático
    # Volume relativo de cadeias alifáticas (A, V, I, L)
    aliphatic_index = (aa_composition['A'] + 2.9 * aa_composition['V']
                      + 3.9 * (aa_composition['I'] + aa_composition['L'])) / length * 100

    # 6. Índice de Instabilidade
    # Prediz estabilidade in vivo baseado em frequência de dipeptídeos
    instability_index = (10 / length) * sum(
        DIWV[sequence[i]][sequence[i+1]] for i in range(length - 1)
    )
    stability = "Estável" if instability_index < 40 else "Instável"

    # 7. Meia-vida (regra do N-terminal)
    # Estimativa grosseira baseada no primeiro aminoácido
    n_term = sequence[0]
    if n_term in ['A','C','G','M','P','S','T','V']:
        half_life = ">10 horas"
    elif n_term in ['I','L','F','W','Y','D','E','N','Q']:
        half_life = "2-30 min"
    else:
        half_life = "Desconhecido"
```

### 3. Formatação de Output

**Formato de texto** (stdout):
```
--- Propriedades Físico-Químicas ---
Comprimento                        : 393
Peso Molecular (Da)                : 43653.11
Ponto Isoelétrico (pI)             : 5.84
GRAVY (Hidropaticidade)            : -0.576
Índice Alifático                   : 78.04
Índice de Instabilidade            : 43.39 (Instável)
Meia-Vida (E. coli, in vitro)      : >10 horas
Total de Resíduos Ácidos (Asp+Glu) : 58
Total de Resíduos Básicos (Arg+Lys+His): 47
Total de Resíduos Polares          : 109
Total de Resíduos Apolares         : 179

--- Composição de Aminoácidos (%) ---
A: 15    (3.82%)   C: 10    (2.54%)   D: 23    (5.85%)
E: 35    (8.91%)   F: 12    (3.05%)   ...
```

**Formato CSV** (`-o output.csv`):
```csv
Propriedade,Valor
Comprimento,393
Peso Molecular (Da),43653.11
Ponto Isoelétrico (pI),5.84
...
```

### 4. Argumentos CLI

**Comando**: `biohub.py physchem`

| Argumento | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `sequence` | Posicional | - | Sequência de aminoácidos |
| `-o, --output` | Opcional | stdout | Salva em CSV |
| `--plot-treemap` | Opcional | - | Treemap de composição (PNG) |
| `--plot-composition` | Opcional | - | Gráfico de barras (PNG) |
| `--plot-hydro` | Opcional | - | Perfil de hidrofobicidade (PNG) |
| `--window` | Opcional | 9 | Janela para perfil hidro |

**Exemplo de uso**:
```bash
python biohub.py physchem MTAMEESQSDISLELPLSQETF \
  -o propriedades.csv \
  --plot-treemap composicao.png \
  --plot-hydro hidrofobicidade.png --window 9
```

### 5. Detalhes de Implementacao de Cada Calculo

#### Argumento: `sequence`

**Parser CLI**: [biohub.py:1025](../biohub.py#L1025)
**Handler**: [biohub.py:371-456](../biohub.py#L371-L456)

**Como funciona**:

1. Conversao e validacao (linhas 373-376):
   - Converte para maiusculo: `sequence = args.sequence.upper()`
   - Verifica se nao esta vazia

2. Calculo de composicao (linha 379):
   - Dictionary comprehension: `{aa: sequence.count(aa) for aa in MOLECULAR_WEIGHT.keys()}`
   - Conta quantas vezes cada aminoacido aparece na sequencia

#### Calculo: Peso Molecular (linha 382)

**Dicionario usado**: MOLECULAR_WEIGHT [biohub.py:38-43](../biohub.py#L38-L43)

**Formula**:
```python
mw = sum(MOLECULAR_WEIGHT[aa] * count for aa, count in aa_composition.items()) - (length - 1) * 18.015
```

**Como funciona**:
1. Soma o peso de todos os aminoacidos
2. Subtrai o peso da agua perdida em cada ligacao peptidica
3. N aminoacidos formam N-1 ligacoes peptidicas
4. Cada ligacao libera 1 H2O (18.015 Da)

**Exemplo**:
- Sequencia: "AAA" (3 alaninas)
- Peso de cada ALA: 89.09 Da
- Peso total: 3 × 89.09 = 267.27 Da
- Ligacoes peptidicas: 2
- Agua perdida: 2 × 18.015 = 36.03 Da
- Peso molecular final: 267.27 - 36.03 = 231.24 Da

#### Calculo: GRAVY (linha 385)

**Dicionario usado**: KYTE_DOOLITTLE [biohub.py:44-49](../biohub.py#L44-L49)

**Formula**:
```python
gravy = sum(KYTE_DOOLITTLE[aa] for aa in sequence) / length
```

**Como funciona**:
1. Para cada aminoacido na sequencia, busca seu valor de hidropaticidade
2. Soma todos os valores
3. Divide pelo comprimento da sequencia (media aritmetica)

**Interpretacao**:
- GRAVY > 0: Proteina hidrofobica (tende a estar no interior ou membrana)
- GRAVY < 0: Proteina hidrofilica (tende a estar exposta ao solvente)
- Range tipico: -2.0 a +2.0

#### Calculo: Ponto Isoeletrico (pI) (linhas 387-400)

**Dicionario usado**: PKA_VALUES [biohub.py:50-54](../biohub.py#L50-L54)

**Algoritmo**:
1. Busca exaustiva de pH 0.00 a 14.00 com passo de 0.01 (1401 iteracoes)
2. Para cada pH, calcula a carga liquida da proteina
3. Guarda o pH onde |carga liquida| eh minima

**Formula de carga liquida** (equacao de Henderson-Hasselbalch):
```python
# Carga do N-terminal (sempre positiva)
net_charge = (10**pKa_N) / (10**pKa_N + 10**pH)

# Carga do C-terminal (sempre negativa)
net_charge -= (10**pH) / (10**pKa_C + 10**pH)

# Para cada aminoacido ionizavel:
# Basicos (R, H, K): contribuem positivamente
net_charge += count * (10**pKa) / (10**pKa + 10**pH)

# Acidos (D, E, C, Y): contribuem negativamente
net_charge -= count * (10**pH) / (10**pKa + 10**pH)
```

**Grupos ionizaveis**:
- N-terminal: pKa = 8.0 (base)
- C-terminal: pKa = 3.65 (acido)
- Asp (D): pKa = 3.9 (acido)
- Glu (E): pKa = 4.07 (acido)
- His (H): pKa = 6.5 (base)
- Cys (C): pKa = 8.5 (acido)
- Tyr (Y): pKa = 10.0 (acido)
- Lys (K): pKa = 10.0 (base)
- Arg (R): pKa = 12.0 (base)

#### Calculo: Indice Alifatico (linha 404)

**Formula**:
```python
aliphatic_index = (A + 2.9*V + 3.9*(I + L)) / length * 100
```

**Como funciona**:
1. Pondera o numero de residuos alifaticos pelo volume relativo
2. A (Ala): peso 1.0
3. V (Val): peso 2.9
4. I (Ile) e L (Leu): peso 3.9
5. Divide pelo comprimento e multiplica por 100

**Interpretacao**:
- Valor alto (>80): Proteina termoestavel
- Relacionado com estabilidade termica
- Proteinas mesofitas: ~50-80
- Proteinas termofitas: >80

#### Calculo: Indice de Instabilidade (linha 407)

**Dicionario usado**: DIWV [biohub.py:60-82](../biohub.py#L60-L82)

**Formula**:
```python
instability_index = (10 / length) * sum(DIWV[seq[i]][seq[i+1]] for i in range(length-1))
```

**Como funciona**:
1. Para cada dipeptideo (par de aminoacidos adjacentes):
   - Busca o peso DIWV[AA1][AA2] na matriz
   - DIWV = Dipeptide Instability Weight Values
2. Soma todos os pesos
3. Multiplica por 10/L (L = comprimento)

**Interpretacao**:
- II < 40: Proteina estavel in vivo
- II >= 40: Proteina instavel in vivo
- Baseado em dados experimentais de E. coli

**Exemplo**:
- Sequencia: "MET-ALA"
- DIWV['M']['A'] = valor da matriz
- Se sequencia tem muitos dipeptideos "raros" ou "desfavoraveis", II aumenta

#### Calculo: Meia-vida (linhas 410-412)

**Regra do N-terminal** (N-end rule):

```python
n_term = sequence[0]
if n_term in ['A','C','G','M','P','S','T','V']:
    half_life = ">10 horas"
elif n_term in ['I','L','F','W','Y','D','E','N','Q']:
    half_life = "2-30 min"
else:
    half_life = "Desconhecido"
```

**Como funciona**:
1. Extrai o primeiro aminoacido da sequencia
2. Consulta tabela empirica baseada no aminoacido N-terminal
3. Estimativa para E. coli in vitro

**Grupos**:
- **Estabilizadores** (>10h): Ala, Cys, Gly, Met, Pro, Ser, Thr, Val
- **Desestabilizadores** (2-30min): Ile, Leu, Phe, Trp, Tyr, Asp, Glu, Asn, Gln
- **Outros**: Arg, Lys, His (variaveis)

**Limitacao**: Estimativa muito grosseira, nao considera:
- Modificacoes pos-traducionais
- Estrutura 3D
- Sequencia interna
- Sistema biologico especifico

#### Flag: `--plot-treemap`

**Parser CLI**: [biohub.py:1027](../biohub.py#L1027)
**Verificacao**: [biohub.py:439-443](../biohub.py#L439-L443)
**Funcao chamada**: `plot_aa_composition_treemap()` em [biohub_viz.py:107-207](../biohub_viz.py#L107-L207)

**Como funciona**:
1. Verifica se modulo biohub_viz esta disponivel (linha 440)
2. Se nao estiver, exibe aviso sobre dependencias (matplotlib, numpy, squarify)
3. Chama funcao de visualizacao passando:
   - aa_composition: dicionario com contagem de cada AA
   - length: comprimento da sequencia
   - output_file: caminho para salvar PNG

#### Flag: `--plot-hydro`

**Parser CLI**: [biohub.py:1029](../biohub.py#L1029)
**Funcao chamada**: `plot_hydropathy_profile()` em [biohub_viz.py:281-362](../biohub_viz.py#L281-L362)

**Como funciona**:
1. Usa janela deslizante de tamanho especificado (--window, padrao 9)
2. Para cada posicao i:
   - Extrai substrings de tamanho window: seq[i:i+window]
   - Calcula media de Kyte-Doolittle da janela
   - Plota score vs posicao
3. Preenche areas:
   - Vermelho: hidrofobico (score > 0)
   - Azul: hidrofilico (score < 0)

**Parametro --window**:
- Padrao: 9 (recomendado por Kyte & Doolittle, 1982)
- Janela maior: perfil mais suave (menos detalhes)
- Janela menor: perfil mais ruidoso (mais detalhes)
- Janela 19-21: ideal para detectar helices transmembrana

---

## Visualizações (biohub_viz.py)

### Treemap de Composição: [biohub_viz.py:107-207](../biohub_viz.py#L107-L207)

```python
def plot_aa_composition_treemap(aa_composition, sequence_length, output_file):
    """Treemap hierárquico agrupado por propriedade química."""

    # Grupos de aminoácidos
    groups = {
        'Hidrofóbico': ['A', 'V', 'L', 'I', 'M', 'F', 'W', 'P'],
        'Polar': ['S', 'T', 'N', 'Q', 'Y', 'C'],
        'Ácido': ['D', 'E'],
        'Básico': ['K', 'R', 'H'],
        'Glicina': ['G']
    }

    # Usa squarify para criar quadrados proporcionais
    squarify.plot(sizes=counts, label=labels, color=colors, ...)
```

### Perfil de Hidrofobicidade: [biohub_viz.py:281-362](../biohub_viz.py#L281-L362)

```python
def plot_hydropathy_profile(sequence, kyte_doolittle_scale, output_file, window=9):
    """Perfil Kyte-Doolittle com janela deslizante."""

    # Calcula média em janela deslizante
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i+window]
        score = sum(kyte_doolittle_scale.get(aa, 0) for aa in window_seq) / window
        scores.append(score)

    # Plota linha com preenchimento (hidrofóbico=vermelho, hidrofílico=azul)
    ax.fill_between(positions, scores, 0,
                     where=[s > 0 for s in scores], color='red')
    ax.fill_between(positions, scores, 0,
                     where=[s < 0 for s in scores], color='blue')
```

---

## 03 - csv2fasta: Conversão CSV → FASTA

### 1. Dicionários e Constantes

Nenhum dicionário específico.

### 2. Função de Cálculo

**Arquivo**: [biohub.py:332-369](../biohub.py#L332-L369)

```python
def handle_csv_to_fasta(args):
    """Converte CSV (ID, Sequência) para FASTA."""

    with open(args.csv_file, mode='r') as infile:
        reader = csv.reader(infile, delimiter=args.delimiter)

        # Se tem cabeçalho, procura colunas por nome ou índice
        if args.header:
            header_row = next(reader)
            try:
                # Tenta converter para inteiro (índice)
                id_col_idx = int(args.id_col)
                seq_col_idx = int(args.seq_col)
            except ValueError:
                # Se não for inteiro, busca pelo nome
                id_col_idx = header_row.index(args.id_col)
                seq_col_idx = header_row.index(args.seq_col)
        else:
            # Sem cabeçalho, deve ser índice numérico
            id_col_idx = int(args.id_col)
            seq_col_idx = int(args.seq_col)

        # Gera registros FASTA
        fasta_records = [
            f">{row[id_col_idx].strip()}\n{row[seq_col_idx].strip()}"
            for row in reader if row
        ]
```

### 3. Formatação de Output

**Input CSV**:
```csv
ID,Sequencia,Organismo
P53_HUMAN,MEEPQSDPSVEPPLSQETFSDLWKLLPEN,Homo sapiens
P53_MOUSE,MTAMEESQSDISLELPLSQETFSGLWKLLP,Mus musculus
```

**Output FASTA**:
```fasta
>P53_HUMAN
MEEPQSDPSVEPPLSQETFSDLWKLLPEN
>P53_MOUSE
MTAMEESQSDISLELPLSQETFSGLWKLLP
```

### 4. Argumentos CLI

**Comando**: `biohub.py csv2fasta`

| Argumento | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `csv_file` | Posicional | - | Arquivo CSV de entrada |
| `-o, --output` | Opcional | stdout | Arquivo FASTA de saída |
| `--id-col` | Opcional | 0 | Coluna do ID (nome ou índice) |
| `--seq-col` | Opcional | 1 | Coluna da sequência (nome ou índice) |
| `--header` | Flag | False | Primeira linha é cabeçalho |
| `--delimiter` | Opcional | `,` | Separador de campos |

**Exemplo de uso**:
```bash
python biohub.py csv2fasta proteinas.csv \
  -o sequencias.fasta \
  --header \
  --id-col ID \
  --seq-col Sequencia
```

### 5. Detalhes de Implementacao

#### Argumento: `csv_file`

**Parser CLI**: [biohub.py:1016](../biohub.py#L1016)
**Handler**: [biohub.py:332-369](../biohub.py#L332-L369)

**Como funciona**:

1. Abertura do arquivo CSV (linha 335-336):
   - Usa modulo csv.reader com delimitador especificado
   - Suporta diferentes encodings (UTF-8 por padrao)

2. Identificacao de colunas (linhas 339-354):

   **Caso 1: CSV com cabecalho** (--header presente):
   - Le primeira linha como header_row (linha 340)
   - Tenta interpretar id_col e seq_col como indices numericos (linha 342)
   - Se falhar (ValueError), interpreta como NOMES de colunas (linha 344-345)
   - Busca indice da coluna pelo nome usando .index()

   **Caso 2: CSV sem cabecalho**:
   - id_col e seq_col DEVEM ser indices numericos (linha 351)
   - Se nao forem, exibe erro (linha 353)

3. Conversao para FASTA (linha 357):
```python
fasta_records = [
    f">{row[id_col_idx].strip()}\n{row[seq_col_idx].strip().replace(' ', '')}"
    for row in reader if row and len(row) > max(id_col_idx, seq_col_idx)
]
```
   - Para cada linha do CSV:
     - Extrai ID da coluna id_col_idx
     - Extrai sequencia da coluna seq_col_idx
     - Remove espacos em branco (.strip())
     - Remove espacos DENTRO da sequencia (.replace(' ', ''))
     - Formata como ">ID\nSEQUENCIA"

#### Flag: `--header`

**Parser CLI**: [biohub.py:1020](../biohub.py#L1020)
**Tipo**: Flag booleana (action="store_true")

**Como funciona**:
- Se presente: args.header = True
- Se ausente: args.header = False
- Quando True, primeira linha eh tratada como cabecalho (linha 339)
- Permite usar nomes de colunas em vez de indices

**Exemplo com header**:
```csv
ID,Sequencia,Organismo
P53_HUMAN,MEEPQSDPSVEPPLSQETFSDLWKLLPEN,Homo sapiens
P53_MOUSE,MTAMEESQSDISLELPLSQETFSGLWKLLP,Mus musculus
```

Comando:
```bash
python biohub.py csv2fasta data.csv --header --id-col ID --seq-col Sequencia
```

**Exemplo sem header**:
```csv
P53_HUMAN,MEEPQSDPSVEPPLSQETFSDLWKLLPEN,Homo sapiens
P53_MOUSE,MTAMEESQSDISLELPLSQETFSGLWKLLP,Mus musculus
```

Comando:
```bash
python biohub.py csv2fasta data.csv --id-col 0 --seq-col 1
```

#### Flag: `--delimiter`

**Parser CLI**: [biohub.py:1021](../biohub.py#L1021)
**Uso**: [biohub.py:336](../biohub.py#L336)
**Padrao**: `,` (virgula)

**Como funciona**:
- Define o caractere separador de campos no CSV
- Passado para csv.reader: `csv.reader(infile, delimiter=args.delimiter)`

**Valores comuns**:
- `,` : CSV padrao (comma-separated values)
- `\t` : TSV (tab-separated values)
- `;` : Comum em locales europeus
- `|` : Pipe-separated

**Exemplo TSV**:
```bash
python biohub.py csv2fasta data.tsv --delimiter $'\t' --header
```

#### Flags: `--id-col` e `--seq-col`

**Parser CLI**: [biohub.py:1018-1019](../biohub.py#L1018-L1019)
**Padroes**: id_col=0, seq_col=1

**Como funciona**:

1. Podem ser especificados como:
   - **Indices numericos** (base 0): `--id-col 0 --seq-col 1`
   - **Nomes de colunas** (requer --header): `--id-col ID --seq-col Sequencia`

2. Logica de interpretacao (linhas 341-345):
```python
try:
    id_col_idx = int(args.id_col)  # Tenta converter para int
except ValueError:
    # Se falhar, busca pelo nome (requer --header)
    id_col_idx = header_row.index(args.id_col)
```

3. Validacao (linhas 346-348):
   - Se coluna nao for encontrada, exibe erro
   - Verifica se cada linha tem colunas suficientes (linha 357)

**Exemplos**:
```bash
# Usando indices (coluna 0 e 1)
python biohub.py csv2fasta data.csv --id-col 0 --seq-col 1

# Usando nomes (requer --header)
python biohub.py csv2fasta data.csv --header --id-col ID --seq-col Sequencia

# CSV com colunas fora de ordem
python biohub.py csv2fasta data.csv --header --id-col Protein --seq-col Seq
```

---

## 04 - contacts: Análise de Contatos Intramoleculares

### 1. Dicionários e Constantes

Nenhum dicionário específico. Usa threshold de distância (padrão: 8.0 Å).

### 2. Função de Cálculo

**Arquivo**: [biohub.py:458-514](../biohub.py#L458-L514)

```python
def calculate_intramolecular_contacts(args):
    """Calcula contatos entre resíduos baseado em distância mínima entre átomos."""

    # Agrupa átomos por resíduo
    residue_atoms = {}  # {res_num: [lista de coordenadas]}

    with open(args.pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                res_num = int(line[22:26])
                coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

                if res_num not in residue_atoms:
                    residue_atoms[res_num] = []
                residue_atoms[res_num].append(coords)

    # Calcula distância mínima entre todos os pares de resíduos
    contacts = []
    for i, r1 in enumerate(residues):
        for r2 in residues[i+1:]:
            # Ignora vizinhos diretos na sequência
            if abs(r1 - r2) <= 1:
                continue

            # Distância mínima entre qualquer átomo de r1 e r2
            min_dist = float('inf')
            for atom1 in residue_atoms[r1]:
                for atom2 in residue_atoms[r2]:
                    # Distância euclidiana 3D
                    dist = sqrt((x1-x2)² + (y1-y2)² + (z1-z2)²)
                    if dist < min_dist:
                        min_dist = dist

            # Se distância < threshold, é um contato
            if min_dist <= args.threshold:
                contacts.append((r1, r2, round(min_dist, 3)))
```

### 3. Formatação de Output

**Formato de texto** (stdout):
```
--- Contatos Intramoleculares (Limiar = 8.0 Å) ---
Res 12 - Res 45: 6.234 Å
Res 12 - Res 89: 7.891 Å
Res 23 - Res 56: 5.678 Å
...
```

**Formato CSV**:
```csv
Residuo_1,Residuo_2,Distancia_A
12,45,6.234
12,89,7.891
23,56,5.678
```

**Visualização - Contact Map**: [biohub_viz.py:545-628](../biohub_viz.py#L545-L628)
```python
def plot_contact_map(contacts, max_residue, output_file, threshold):
    """Mapa de contatos como heatmap (matriz de distâncias)."""

    # Cria matriz NxN (N = número de resíduos)
    contact_matrix = np.zeros((max_residue, max_residue))

    # Preenche matriz com distâncias
    for res1, res2, dist in contacts:
        contact_matrix[res1-1, res2-1] = dist
        contact_matrix[res2-1, res1-1] = dist  # Simétrica

    # Heatmap com escala de cores (viridis_r)
    ax.imshow(contact_matrix, cmap='viridis_r', ...)
```

### 4. Argumentos CLI

**Comando**: `biohub.py contacts`

| Argumento | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `pdb_file` | Posicional | - | Arquivo PDB de entrada |
| `-t, --threshold` | Opcional | 8.0 | Distância máxima (Å) |
| `-o, --output` | Opcional | stdout | Salva em CSV |
| `--plot` | Opcional | - | Mapa de contatos (PNG) |

**Exemplo de uso**:
```bash
python biohub.py contacts 1TUP.pdb \
  -t 8.0 \
  -o contatos.csv \
  --plot contact_map.png
```

### 5. Detalhes de Implementacao

#### Argumento: `pdb_file`

**Parser CLI**: [biohub.py:1034](../biohub.py#L1034)
**Handler**: [biohub.py:458-514](../biohub.py#L458-L514)

**Como funciona**:

1. Leitura e agrupamento de atomos (linhas 465-476):
   - Le arquivo PDB linha por linha
   - Para cada linha ATOM da primeira cadeia:
     - Extrai res_num (colunas 22-26)
     - Extrai coordenadas x, y, z (colunas 30-38, 38-46, 46-54)
     - Agrupa todos os atomos do mesmo residuo em uma lista
   - Estrutura: `residue_atoms = {res_num: [(x,y,z), (x,y,z), ...]}`

2. Calculo de distancias minimas (linhas 482-498):
   - Para cada par de residuos (r1, r2):
     - Ignora pares adjacentes na sequencia (|r1-r2| <= 1) (linha 485)
     - Calcula distancia entre TODOS os pares de atomos
     - Guarda a distancia MINIMA entre qualquer atomo de r1 e r2

3. Formula da distancia euclidiana (linha 492):
```python
dist = math.sqrt(sum((c1-c2)**2 for c1,c2 in zip(atom1, atom2)))
# Equivalente a: sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
```

**Por que distancia minima?**
- Residuos podem estar em contato mesmo se C-alfa estao longe
- Exemplo: cadeias laterais longas (Lys, Arg) podem interagir
- Considera geometria real da interacao

**Por que ignora vizinhos diretos?**
- Residuos adjacentes (i, i+1) sempre estao proximos por definicao
- Ligacao peptidica conecta diretamente
- Foco em contatos nao-locais (estrutura terciaria)

#### Flag: `-t, --threshold`

**Parser CLI**: [biohub.py:1035](../biohub.py#L1035)
**Uso**: [biohub.py:497](../biohub.py#L497)
**Padrao**: 8.0 Angstrom

**Como funciona**:

1. Compara distancia minima com threshold (linha 497):
```python
if min_dist <= args.threshold:
    contacts.append((r1, r2, round(min_dist, 3)))
```

2. Apenas pares com distancia <= threshold sao considerados contatos

**Valores tipicos**:
- 4.0 A: Contatos muito fortes (ponte de hidrogenio, ponte salina)
- 6.0 A: Contatos fortes (interacoes Van der Waals)
- 8.0 A: Threshold padrao (cobre maioria das interacoes nao-covalentes)
- 10.0 A: Contatos fracos ou proximos

**Interpretacao**:
- Threshold baixo: menos contatos, mais especificos
- Threshold alto: mais contatos, menos especificos

#### Flag: `--plot`

**Parser CLI**: [biohub.py:1037](../biohub.py#L1037)
**Verificacao**: [biohub.py:509-514](../biohub.py#L509-L514)
**Funcao chamada**: `plot_contact_map()` em [biohub_viz.py:545-628](../biohub_viz.py#L545-L628)

**Como funciona**:

1. Gera mapa de contatos (contact map) como heatmap
2. Cria matriz NxN onde N = numero de residuos
3. Preenche matriz com distancias entre residuos
4. Matriz eh simetrica: dist(i,j) = dist(j,i)
5. Visualizacao:
   - Eixo X: residuo i
   - Eixo Y: residuo j
   - Cor: distancia (escala viridis invertida)
   - Pontos escuros: contatos proximos
   - Pontos claros: residuos distantes

---

## 05 - hydrophoby: Perfil de Hidrofobicidade por Resíduo

### 1. Dicionários e Constantes

**KYTE_DOOLITTLE**: [biohub.py:44-49](../biohub.py#L44-L49) (mesmo usado em physchem)
**THREE_TO_ONE**: [biohub.py:31-36](../biohub.py#L31-L36)

### 2. Função de Cálculo

**Arquivo**: [biohub.py:719-788](../biohub.py#L719-L788)

```python
def predict_solvent_hydrophoby(args):
    """Atribui hidrofobicidade Kyte-Doolittle a cada átomo."""

    atoms = parse_pdb_atoms(args.pdb_file)
    results = []

    # Para cada átomo, atribui hidrofobicidade do resíduo
    for atom in atoms:
        res_name = atom["res_name"]  # Ex: "ALA"

        # Converte 3 letras → 1 letra
        aa_code = THREE_TO_ONE.get(res_name, None)  # Ex: "A"

        if aa_code:
            # Busca score na escala Kyte-Doolittle
            hydro_score = KYTE_DOOLITTLE.get(aa_code, 0.0)
        else:
            # Se não for aminoácido padrão (ligante), define 0
            hydro_score = 0.0

        results.append({
            "chain_id": atom["chain_id"],
            "res_num": atom["res_num"],
            "atom_num": atom["atom_num"],
            "hydrophobicity": hydro_score
        })
```

**Função auxiliar - Parser de PDB**: [biohub.py:111-134](../biohub.py#L111-L134)
```python
def parse_pdb_atoms(pdb_filepath: str):
    """Extrai coordenadas e informações de todos os átomos."""
    atoms = []
    with open(pdb_filepath, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append({
                    "atom_num": int(line[6:11]),
                    "atom_name": line[12:16].strip(),
                    "res_name": line[17:20].strip(),
                    "chain_id": line[21],
                    "res_num": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": line[76:78].strip().upper()
                })
    return atoms
```

### 3. Formatação de Output

**Formato CSV**:
```csv
Chain,ResNum,ResName,AtomNum,AtomName,Hydrophobicity
A,12,ALA,123,CA,1.800
A,12,ALA,124,CB,1.800
A,13,ARG,125,N,-4.500
```

**PDB Anotado** (`--write-pdb`): [biohub.py:516-564](../biohub.py#L516-L564)
```python
def write_pdb_with_bfactor(input_pdb_path, output_pdb_path, atom_values, property_name):
    """Reescreve PDB substituindo B-factor por valor calculado."""

    with open(input_pdb_path, 'r') as infile, open(output_pdb_path, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_num = int(line[6:11])

                # Busca valor calculado para este átomo
                new_bfactor = value_by_atom[atom_num]

                # Substitui colunas 61-66 (B-factor)
                new_line = line[:60] + f"{new_bfactor:6.2f}" + line[66:]
                outfile.write(new_line)
```

**Sessão PyMOL** (`--pymol`): [biohub.py:566-717](../biohub.py#L566-L717)
```python
def generate_pymol_session(pdb_path, output_pse, property_type="hydrophobicity", ...):
    """Gera script .pml e sessão .pse com visualização configurada."""

    # Script PyMOL gerado:
    """
    load proteina.pdb
    show cartoon
    show surface
    set transparency, 0.5
    spectrum b, blue_white_red, minimum=-4.5, maximum=4.5
    """
    # Azul = hidrofílico (-4.5), Vermelho = hidrofóbico (+4.5)
```

**Visualização - Perfil**: [biohub_viz.py:365-458](../biohub_viz.py#L365-L458)
```python
def plot_hydrophoby_profile(hydrophoby_per_atom, output_file):
    """Gráfico de linha mostrando hidrofobicidade por posição."""

    # Calcula média por resíduo (todos os átomos têm mesmo valor)
    residue_hydro = {}
    for atom in hydrophoby_per_atom:
        residue_key = (atom['chain_id'], atom['res_num'])
        residue_hydro[residue_key] = atom['hydrophobicity']

    # Plot com cores divergentes (azul-branco-vermelho)
    scatter = ax.scatter(residue_numbers, hydro_values, c=hydro_values,
                         cmap='RdBu_r', vmin=-4.5, vmax=4.5)

    # Preenche áreas
    ax.fill_between(positions, values, 0,
                     where=[h > 0], color='red', label='Hidrofóbico')
    ax.fill_between(positions, values, 0,
                     where=[h < 0], color='blue', label='Hidrofílico')
```

### 4. Argumentos CLI

**Comando**: `biohub.py hydrophoby`

| Argumento | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `pdb_file` | Posicional | - | Arquivo PDB de entrada |
| `-o, --output` | Opcional | stdout | Salva em CSV |
| `--write-pdb` | Opcional | - | PDB com hidro no B-factor |
| `--pymol` | Opcional | - | Sessão PyMOL (.pse) |
| `--plot-hydrophoby` | Opcional | - | Gráfico de perfil (PNG) |

**Exemplo de uso**:
```bash
python biohub.py hydrophoby 1TUP.pdb \
  -o hidrofobicidade.csv \
  --write-pdb 1TUP_hydro.pdb \
  --pymol visualizacao_hydro.pse \
  --plot-hydrophoby perfil_hydro.png
```

### 5. Detalhes de Implementacao

#### Argumento: `pdb_file`

**Parser CLI**: [biohub.py:1041](../biohub.py#L1041)
**Handler**: [biohub.py:719-788](../biohub.py#L719-L788)
**Funcao auxiliar**: `parse_pdb_atoms()` em [biohub.py:111-134](../biohub.py#L111-L134)

**Como funciona**:

1. Parse do arquivo PDB (linhas 721-722):
   - Chama `parse_pdb_atoms()` que le todas as linhas ATOM e HETATM
   - Extrai para cada atomo:
     - atom_num (colunas 6-11)
     - atom_name (colunas 12-16)
     - res_name (colunas 17-20)
     - chain_id (coluna 21)
     - res_num (colunas 22-26)
     - coordenadas x, y, z (colunas 30-38, 38-46, 46-54)
     - element (colunas 76-78)

2. Atribuicao de hidrofobicidade (linhas 727-746):
   - Para cada atomo:
     - Extrai res_name (ex: "ALA", "VAL")
     - Converte para codigo 1 letra usando THREE_TO_ONE (ex: "A", "V")
     - Busca score no dicionario KYTE_DOOLITTLE
     - TODOS os atomos do mesmo residuo recebem o MESMO score
     - Se nao for aminoacido padrao (ligante), score = 0.0

**Por que atribuir por residuo?**
- Hidrofobicidade eh propriedade do aminoacido, nao do atomo
- Facilita visualizacao em PyMOL/VMD
- Permite colorir toda a cadeia lateral uniformemente

#### Flag: `--write-pdb`

**Parser CLI**: [biohub.py:1043](../biohub.py#L1043)
**Uso**: [biohub.py:774-777](../biohub.py#L774-L777)
**Funcao chamada**: `write_pdb_with_bfactor()` em [biohub.py:516-564](../biohub.py#L516-L564)

**Como funciona**:

1. Prepara dados (linha 776):
```python
atom_values = [{'atom_num': atom['atom_num'], 'value': atom['hydrophobicity']} for atom in results]
```

2. Funcao `write_pdb_with_bfactor()`:
   - Le arquivo PDB original linha por linha
   - Para linhas ATOM/HETATM:
     - Extrai atom_num (colunas 6-11)
     - Busca valor calculado no dicionario atom_values
     - Substitui colunas 61-66 (B-factor) pelo novo valor
     - Formato: `f"{valor:6.2f}"` (6 caracteres, 2 decimais)
   - Mantem todas as outras linhas inalteradas

**Formato B-factor no PDB**:
```
ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N
                                                           ^^^^^^
                                                       Colunas 61-66
                                                       (B-factor original)
```

Apos processamento:
```
ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00  1.90           N
                                                           ^^^^^^
                                                        Hidrofobicidade
```

#### Flag: `--pymol`

**Parser CLI**: [biohub.py:1044](../biohub.py#L1044)
**Uso**: [biohub.py:780-781](../biohub.py#L780-L781)
**Funcao chamada**: `generate_pymol_session()` em [biohub.py:566-717](../biohub.py#L566-L717)

**Como funciona**:

1. Verifica dependencias (linha 780):
   - Requer --write-pdb (arquivo PDB com valores no B-factor)
   - Se --pymol sem --write-pdb, cria PDB temporario

2. Parametros passados (linha 781):
   - property_type="hydrophobicity"
   - min_val=-4.5 (Arg, mais hidrofilico)
   - max_val=4.5 (Ile, mais hidrofobico)

3. Script PyMOL gerado:
```python
load arquivo.pdb
show cartoon
show surface
set transparency, 0.5
spectrum b, blue_white_red, minimum=-4.5, maximum=4.5
```

**Esquema de cores**:
- Azul: hidrofilico (B-factor = -4.5)
- Branco: neutro (B-factor = 0.0)
- Vermelho: hidrofobico (B-factor = +4.5)

4. Salva sessao .pse:
   - Executa PyMOL em modo batch
   - Aplica comandos do script
   - Salva estado como arquivo .pse
   - Usuario pode abrir diretamente no PyMOL

---

## 06 - sasa: Área de Superfície Acessível ao Solvente

### 1. Dicionários e Constantes

**VDW_RADII** (Raios de Van der Waals): [biohub.py:55-59](../biohub.py#L55-L59)
```python
VDW_RADII = {
    'H': 1.20,   # Hidrogênio
    'C': 1.70,   # Carbono
    'N': 1.55,   # Nitrogênio
    'O': 1.52,   # Oxigênio
    'S': 1.80,   # Enxofre
    'P': 1.80,   # Fósforo
    'DEFAULT': 1.70  # Valor padrão
}
```

**Constantes de cálculo**:
- **Probe radius** (raio da sonda de água): 1.4 Å
- **Número de pontos na esfera** (Fibonacci sphere): 960 (padrão)

### 2. Função de Cálculo

**Arquivo**: [biohub.py:790-943](../biohub.py#L790-L943)

#### Algoritmo: Shrake-Rupley

```python
def calculate_sasa(args):
    """Calcula SASA usando método de Shrake-Rupley."""

    atoms = parse_pdb_atoms(args.pdb_file)

    # 1. Gera pontos uniformemente distribuídos na esfera (Fibonacci sphere)
    sphere_points = generate_sphere_points(args.num_points)  # Ex: 960 pontos

    total_sasa = 0.0
    sasa_per_atom = []

    # 2. Para cada átomo...
    for i, atom_i in enumerate(atoms):
        # Raio estendido = raio VdW + raio da sonda (água = 1.4Å)
        radius_i = VDW_RADII.get(atom_i["element"], 1.70)
        extended_radius = radius_i + args.probe_radius

        accessible_points = 0

        # 3. Testa cada ponto na superfície estendida
        for sp in sphere_points:
            # Coordenadas do ponto de teste
            point = (
                atom_i["x"] + extended_radius * sp[0],
                atom_i["y"] + extended_radius * sp[1],
                atom_i["z"] + extended_radius * sp[2]
            )

            point_is_accessible = True

            # 4. Verifica se ponto está ocluído por outro átomo
            for j, atom_j in enumerate(atoms):
                if i == j:
                    continue

                radius_j = VDW_RADII.get(atom_j["element"], 1.70)

                # Distância do ponto ao centro do átomo j
                dist_squared = sum((p - c)**2 for p, c in zip(
                    point, (atom_j["x"], atom_j["y"], atom_j["z"])
                ))

                # Se ponto está dentro da esfera estendida de j, está ocluído
                if dist_squared < (radius_j + args.probe_radius)**2:
                    point_is_accessible = False
                    break

            if point_is_accessible:
                accessible_points += 1

        # 5. SASA = proporção de pontos acessíveis × área da esfera
        atom_sasa = (accessible_points / args.num_points) * 4π * extended_radius²
        total_sasa += atom_sasa

        sasa_per_atom.append({
            "atom_num": atom_i["atom_num"],
            "res_num": atom_i["res_num"],
            "sasa": atom_sasa
        })
```

**Função auxiliar - Fibonacci Sphere**: [biohub.py:197-208](../biohub.py#L197-L208)
```python
def generate_sphere_points(n_points: int):
    """Distribui n_points uniformemente na superfície de esfera unitária."""
    points = []
    phi = (1 + sqrt(5)) / 2  # Proporção áurea

    for i in range(n_points):
        y = 1 - (2 * i / (n_points - 1))
        radius = sqrt(1 - y²)
        theta = 2π * i / phi

        points.append((
            cos(theta) * radius,  # x
            y,                     # y
            sin(theta) * radius   # z
        ))
    return points
```

### 3. Formatação de Output

**Formato CSV**:
```csv
Chain,ResNum,ResName,AtomNum,AtomName,SASA_A2
A,12,ALA,123,N,0.00
A,12,ALA,124,CA,2.34
A,12,ALA,125,C,0.00
A,12,ALA,126,O,15.67
```

**Terminal**:
```
SASA Total da Molécula: 8542.34 Å²
Total de átomos analisados: 3045

--- SASA por Átomo ---
Chain | ResNum | ResName | AtomNum | AtomName | SASA (Ų)
A     | 12     | ALA     | 123     | N        | 0.00
A     | 12     | ALA     | 124     | CA       | 2.34
...
```

**PDB Anotado** (`--write-pdb`):
- Calcula **SASA médio por resíduo** (mais relevante biologicamente)
- Escreve média no B-factor de todos os átomos do resíduo
- Usa percentil 70 para range de visualização

**Visualização - Perfil**: [biohub_viz.py:461-542](../biohub_viz.py#L461-L542)
```python
def plot_sasa_profile(sasa_per_atom, output_file):
    """Gráfico de linha mostrando SASA médio por resíduo."""

    # Calcula média por resíduo
    residue_sasa = defaultdict(list)
    for atom in sasa_per_atom:
        residue_key = (atom['chain_id'], atom['res_num'])
        residue_sasa[residue_key].append(atom['sasa'])

    avg_sasa = [sum(sasa_list)/len(sasa_list) for sasa_list in residue_sasa.values()]

    # Threshold típico: 20 Ų
    threshold = 20.0

    # Plot com cores (exposto=azul, enterrado=vermelho)
    scatter = ax.scatter(residue_numbers, avg_sasa, c=avg_sasa,
                         cmap='RdYlBu_r')  # Vermelho-Amarelo-Azul invertido

    # Linha de threshold
    ax.axhline(y=threshold, color='red', linestyle='--',
               label=f'Threshold ({threshold:.0f} Ų)')

    # Preenche áreas
    ax.fill_between(positions, avg_sasa, threshold,
                     where=[s > threshold], color='blue', label='Exposto')
    ax.fill_between(positions, avg_sasa, threshold,
                     where=[s <= threshold], color='red', label='Enterrado')
```

### 4. Argumentos CLI

**Comando**: `biohub.py sasa`

| Argumento | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `pdb_file` | Posicional | - | Arquivo PDB de entrada |
| `--probe-radius` | Opcional | 1.4 | Raio da sonda (Å) |
| `--num-points` | Opcional | 960 | Pontos na esfera Fibonacci |
| `-o, --output` | Opcional | stdout | Salva em CSV |
| `--write-pdb` | Opcional | - | PDB com SASA no B-factor |
| `--pymol` | Opcional | - | Sessão PyMOL (.pse) |
| `--plot-profile` | Opcional | - | Gráfico de perfil (PNG) |

**Exemplo de uso**:
```bash
python biohub.py sasa 1TUP.pdb \
  --probe-radius 1.4 \
  --num-points 960 \
  -o sasa.csv \
  --write-pdb 1TUP_sasa.pdb \
  --pymol visualizacao_sasa.pse \
  --plot-profile perfil_sasa.png
```

**Interpretação**:
- **SASA > 20 Ų**: Resíduo exposto ao solvente
- **SASA ≤ 20 Ų**: Resíduo enterrado no interior
- **SASA total**: Indica compactação da estrutura

### 5. Detalhes de Implementacao - Algoritmo de Shrake-Rupley

#### Argumento: `pdb_file`

**Parser CLI**: [biohub.py:1049](../biohub.py#L1049)
**Handler**: [biohub.py:790-943](../biohub.py#L790-L943)
**Algoritmo**: Shrake-Rupley (1973)

**Como funciona - Etapa por etapa**:

1. Parse de atomos (linha 792):
   - Usa `parse_pdb_atoms()` para extrair todos os atomos
   - Obtem coordenadas x,y,z e elemento quimico de cada atomo

2. Geracao de pontos na esfera (linha 795):
   - Chama `generate_sphere_points(num_points)` [biohub.py:197-208](../biohub.py#L197-L208)
   - Algoritmo: Fibonacci Sphere
   - Distribui uniformemente num_points na superficie de esfera unitaria
   - Pontos servem como "direcoes" para testar acessibilidade

**Algoritmo Fibonacci Sphere**:
```python
phi = (1 + sqrt(5)) / 2  # Proporcao aurea
for i in range(n_points):
    y = 1 - (2 * i / (n_points - 1))
    radius = sqrt(1 - y**2)
    theta = 2*pi * i / phi
    points.append((cos(theta)*radius, y, sin(theta)*radius))
```

3. Para CADA atomo i (linhas 800-829):

   **Passo 3.1** - Calcula raio estendido (linhas 801-802):
   ```python
   radius_i = VDW_RADII[elemento]  # Ex: C=1.70, N=1.55, O=1.52
   extended_radius = radius_i + probe_radius  # Ex: 1.70 + 1.4 = 3.10 A
   ```

   **Passo 3.2** - Para cada ponto na esfera (linhas 805-815):
   - Calcula coordenadas do ponto de teste (linha 807):
     ```python
     point = (atom_x + extended_radius * sp_x,
              atom_y + extended_radius * sp_y,
              atom_z + extended_radius * sp_z)
     ```
   - Verifica se ponto esta ocluido por outros atomos (linhas 809-814):
     - Para cada OUTRO atomo j:
       - Calcula distancia do ponto ao centro do atomo j
       - Se distancia < (radius_j + probe_radius), ponto esta OCLUIDO
       - Se ocluido, break (nao precisa testar outros atomos)
   - Se nao ocluido, ponto esta ACESSIVEL (linha 815)

   **Passo 3.3** - Calcula SASA do atomo (linha 818):
   ```python
   atom_sasa = (accessible_points / num_points) * 4*pi * extended_radius^2
   ```
   - Formula: (fracao acessivel) × (area da esfera)
   - Area da esfera: 4πr²
   - Se todos os pontos acessiveis: SASA = area total da esfera
   - Se metade acessivel: SASA = metade da area

**Exemplo numerico**:
- Atomo: Carbono (radius = 1.70 A)
- Probe radius: 1.4 A
- Extended radius: 3.10 A
- Num points: 960
- Accessible points: 480 (metade)
- Area total: 4π × 3.10² = 120.76 A²
- SASA: (480/960) × 120.76 = 60.38 A²

#### Flag: `--probe-radius`

**Parser CLI**: [biohub.py:1050](../biohub.py#L1050)
**Uso**: [biohub.py:802](../biohub.py#L802)
**Padrao**: 1.4 Angstrom

**Como funciona**:
- Raio da molecula de agua: 1.4 A
- Define o tamanho da "sonda" que testa acessibilidade
- Extended radius = VDW radius + probe radius
- Simula a menor distancia que agua pode se aproximar

**Valores alternativos**:
- 1.2 A: Sonda menor (aumenta SASA calculado)
- 1.4 A: Padrao (agua)
- 1.6 A: Sonda maior (diminui SASA calculado)

#### Flag: `--num-points`

**Parser CLI**: [biohub.py:1051](../biohub.py#L1051)
**Uso**: [biohub.py:795, 818](../biohub.py#L795)
**Padrao**: 960

**Como funciona**:
- Numero de pontos distribuidos na esfera de cada atomo
- Mais pontos = maior precisao, mas mais lento
- Complexidade: O(N² × M) onde N=atomos, M=pontos

**Valores tipicos**:
- 92 pontos: Rapido, baixa precisao
- 240 pontos: Balanceado
- 960 pontos: Padrao recomendado (precisao boa)
- 3840 pontos: Alta precisao, muito lento

**Trade-off**:
- 960 -> 3840: melhora ~2% na precisao, 4x mais lento
- 240 -> 960: melhora ~5% na precisao, 4x mais lento

#### Flag: `--write-pdb`

**Parser CLI**: [biohub.py:1053](../biohub.py#L1053)
**Uso**: [biohub.py:858-892](../biohub.py#L858-L892)

**Como funciona**:

1. Calcula SASA medio por residuo (linhas 860-875):
   - Agrupa todos os atomos de cada residuo
   - Calcula media dos valores SASA
   - Por que media? SASA total por residuo seria muito alto

2. Exemplo:
   - Residuo ALA tem 5 atomos: N, CA, C, O, CB
   - SASA dos atomos: 10, 5, 8, 20, 15 A²
   - SASA medio: (10+5+8+20+15)/5 = 11.6 A²

3. Escreve no B-factor (similar a hydrophoby):
   - TODOS os atomos do mesmo residuo recebem o mesmo valor (a media)
   - Facilita visualizacao uniforme do residuo

#### Flag: `--plot-profile`

**Parser CLI**: [biohub.py:1055](../biohub.py#L1055)
**Funcao chamada**: `plot_sasa_profile()` em [biohub_viz.py:461-542](../biohub_viz.py#L461-L542)

**Como funciona**:
1. Calcula SASA medio por residuo
2. Plota grafico: posicao do residuo (eixo X) vs SASA (eixo Y)
3. Linha horizontal no threshold (20 A²)
4. Preenche areas:
   - Azul: SASA > 20 (exposto ao solvente)
   - Vermelho: SASA <= 20 (enterrado)
5. Colorbar mostra escala de exposicao

**Interpretacao biologica**:
- Residuos expostos (SASA alto):
  - Candidatos para mutagenese dirigida
  - Epitopos para anticorpos
  - Sites de modificacao pos-traducional
- Residuos enterrados (SASA baixo):
  - Nucleo hidrofobico
  - Sites cataliticos protegidos
  - Residuos estruturais

---

## Módulo de Visualização (biohub_viz.py)

### Constantes de Padronização Visual: [biohub_viz.py:38-78](../biohub_viz.py#L38-L78)

```python
# Posicionamento padronizado
LEGEND_X_POSITION = 1.05         # Posição X para legendas à direita
LEGEND_X_WITH_COLORBAR = 1.18    # Para gráficos com colorbar
LEGEND_SPACING = 0.03            # Espaçamento vertical

# Posições Y
LEGEND_TOP_Y = 0.98              # Legenda no topo
INFO_BOX_1_Y = 0.65              # Primeira caixa de info
INFO_BOX_2_Y = 0.40              # Segunda caixa
INFO_BOX_3_Y = 0.15              # Terceira caixa

# Estilo de legendas
LEGEND_STYLE = {
    'frameon': True,
    'shadow': True,
    'fancybox': True,
    'fontsize': 10
}

# Estilo de caixas de informação
INFO_BOX_STYLE = {
    'boxstyle': 'round',
    'facecolor': 'wheat',
    'alpha': 0.9,
    'edgecolor': 'black',
    'linewidth': 1.0
}

# Margens para layout
LAYOUT_MARGINS_WITH_LEGEND = [0, 0, 0.85, 1]     # Com legenda à direita
LAYOUT_MARGINS_WITH_COLORBAR = [0, 0, 0.78, 1]   # Com colorbar
```

### Dependências e Verificação: [biohub_viz.py:12-104](../biohub_viz.py#L12-L104)

```python
# Importações opcionais
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    matplotlib.use('Agg')  # Backend não-interativo
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import squarify
    HAS_SQUARIFY = True
except ImportError:
    HAS_SQUARIFY = False

def check_dependencies(require_numpy=False, require_squarify=False):
    """Verifica se dependências necessárias estão disponíveis."""
    if not HAS_MATPLOTLIB:
        print("Erro: matplotlib não instalado.")
        return False

    if require_numpy and not HAS_NUMPY:
        print("Erro: numpy não instalado.")
        return False

    if require_squarify and not HAS_SQUARIFY:
        print("Aviso: squarify não instalado. Usando alternativa.")

    return True
```

---

## Resumo das Funções por Módulo

### biohub.py - Funções Principais

| Função | Linha | Propósito |
|--------|-------|-----------|
| `print_banner()` | 86-103 | Exibe logo ASCII |
| `parse_pdb_atoms()` | 111-134 | Parser de coordenadas PDB |
| `get_sequence_from_pdb()` | 136-160 | Extrai sequência FASTA |
| `extract_pdb_header_info()` | 162-184 | Extrai metadados do cabeçalho |
| `write_csv()` | 186-195 | Salva dados em CSV |
| `generate_sphere_points()` | 197-208 | Fibonacci sphere para SASA |
| `filter_pdb_content()` | 212-267 | Filtra cadeias/água/ligantes |
| `handle_fetch_pdb()` | 269-316 | Download de PDB |
| `handle_pdb_to_fasta()` | 318-330 | Conversão PDB→FASTA |
| `handle_csv_to_fasta()` | 332-369 | Conversão CSV→FASTA |
| `calculate_physicochemical_properties()` | 371-456 | Propriedades físico-químicas |
| `calculate_intramolecular_contacts()` | 458-514 | Contatos entre resíduos |
| `write_pdb_with_bfactor()` | 516-564 | Anota PDB com valores no B-factor |
| `generate_pymol_session()` | 566-717 | Gera script/sessão PyMOL |
| `predict_solvent_hydrophoby()` | 719-788 | Hidrofobicidade por átomo |
| `calculate_sasa()` | 790-943 | SASA por átomo (Shrake-Rupley) |
| `run_apbs_analysis()` | 945-983 | Energia de solvatação (PDB2PQR+APBS) |
| `main()` | 987-1090 | Parser CLI e dispatcher |

### biohub_viz.py - Funções de Visualização

| Função | Linha | Propósito |
|--------|-------|-----------|
| `check_dependencies()` | 81-104 | Verifica matplotlib/numpy/squarify |
| `plot_aa_composition_treemap()` | 107-207 | Treemap de composição de AA |
| `plot_aa_composition_bar()` | 210-278 | Gráfico de barras de composição |
| `plot_hydropathy_profile()` | 281-362 | Perfil Kyte-Doolittle (janela deslizante) |
| `plot_hydrophoby_profile()` | 365-458 | Perfil de hidrofobicidade por resíduo |
| `plot_sasa_profile()` | 461-542 | Perfil de SASA por resíduo |
| `plot_contact_map()` | 545-628 | Mapa de contatos (heatmap) |

---

## Detalhes Científicos para Apresentação

### 1. Escala de Kyte-Doolittle
- **Range**: -4.5 (mais hidrofílico: Arg) a +4.5 (mais hidrofóbico: Ile)
- **Interpretação**: Valores positivos = resíduos tendem a estar no interior da proteína
- **Janela deslizante**: Tamanho típico = 9 resíduos (média móvel)
- **Threshold transmembrana**: Picos > 1.6 sugerem domínios transmembrana

### 2. Ponto Isoelétrico (pI)
- **Definição**: pH onde carga líquida = 0
- **Método**: Busca por pH entre 0-14 com passos de 0.01
- **Grupos ionizáveis**: N-term, C-term, D, E, H, C, Y, K, R
- **Equação de Henderson-Hasselbalch**: Calcula carga de cada grupo em cada pH

### 3. Índice de Instabilidade
- **Threshold**: < 40 = estável, ≥ 40 = instável
- **Método**: Frequência de dipeptídeos instáveis (matriz DIWV)
- **Fórmula**: II = (10/L) × Σ DIWV[AAᵢ][AAᵢ₊₁]
- **Limitação**: Desenvolvido para E. coli in vitro

### 4. SASA (Shrake-Rupley)
- **Probe radius**: 1.4 Å (raio da molécula de água)
- **Pontos na esfera**: 960 (compromisso precisão/performance)
- **Threshold**: 20 Ų = limiar típico para exposto vs enterrado
- **Complexidade**: O(N² × M) onde N=átomos, M=pontos

### 5. Contatos Intramoleculares
- **Threshold**: 8.0 Å (distância típica para interações não-covalentes)
- **Método**: Distância mínima entre qualquer átomo de dois resíduos
- **Ignora**: Vizinhos diretos na sequência (i, i±1)
- **Aplicação**: Predição de estrutura secundária/terciária

---

## Sugestões para Apresentação

### Slide 1: Visão Geral
- Logo do BioHub
- 6 módulos principais
- Fluxo típico: fetchpdb → fasta → physchem → sasa/hydrophoby/contacts

### Slide 2-7: Um slide por função
Para cada função, mostrar:
1. **Dicionário/Constante** (canto superior esquerdo)
2. **Código do cálculo** (centro, destacado)
3. **Output visual** (direita: gráfico ou tabela)
4. **Exemplo CLI** (rodapé)

### Slide 8: Integração com ferramentas externas
- PyMOL (visualização 3D)
- PDB2PQR + APBS (energia eletrostática)
- Matplotlib/NumPy (gráficos científicos)

### Slide 9: Casos de uso
- Análise rápida de proteína do PDB
- Predição de regiões transmembrana
- Comparação de estabilidade entre mutantes
- Identificação de resíduos de superfície para mutagênese

---

## Referências Científicas

1. **Kyte-Doolittle Scale**: Kyte, J., & Doolittle, R. F. (1982). Journal of Molecular Biology, 157(1), 105-132.

2. **Shrake-Rupley Algorithm**: Shrake, A., & Rupley, J. A. (1973). Journal of Molecular Biology, 79(2), 351-371.

3. **Instability Index**: Guruprasad, K., et al. (1990). Protein Engineering, 4(2), 155-161.

4. **Aliphatic Index**: Ikai, A. (1980). Journal of Biochemistry, 88(6), 1895-1898.

5. **Fibonacci Sphere**: Swinbank, R., & Purser, R. J. (2006). Quarterly Journal of the Royal Meteorological Society, 132(619), 1769-1793.

---

## Conclusão

O BioHub demonstra como:
- **Hardcoded implementations** podem ser educacionais e funcionais
- **Zero dependências externas** para análises básicas (apenas stdlib)
- **Visualizações opcionais** adicionam valor sem comprometer a funcionalidade core
- **CLI bem estruturado** facilita automação e integração com pipelines

**Próximos passos sugeridos**:
- Suporte a estruturas multi-cadeia
- Cálculo de estrutura secundária (DSSP-like)
- Análise de dinâmica molecular (trajetórias)
- Interface web (Flask/Streamlit)


