# BioHub - Guia Detalhado de Funções para Apresentação

Este documento organiza todas as funções do BioHub ([biohub.py](../biohub.py) e [biohub_viz.py](../biohub_viz.py)) de forma didática para apresentação, separando:
1. **Dicionários e Constantes** - Dados científicos necessários
2. **Função de Cálculo** - Implementação do algoritmo
3. **Formatação de Output** - Como os resultados são apresentados
4. **Argumentos CLI** - Parâmetros que o usuário pode usar

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


