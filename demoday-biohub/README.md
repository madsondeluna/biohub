# BioHub - Roteiro Completo de Demonstração

Este roteiro demonstra **TODAS** as funcionalidades do BioHub de ponta a ponta usando o PDB **1TUP** (p53 Tumor Suppressor).

## Índice
1. [Download e Limpeza do PDB](#1-download-e-limpeza-do-pdb)
2. [Conversão PDB para FASTA](#2-conversão-pdb-para-fasta)
3. [Análises com FASTA](#3-análises-com-fasta)
4. [Conversão CSV para FASTA](#4-conversão-csv-para-fasta)
5. [Análise de Contatos Intramoleculares](#5-análise-de-contatos-intramoleculares)
6. [Análise de Exposição ao Solvente (Predição)](#6-análise-de-exposição-ao-solvente-predição)
7. [Cálculo de SASA (Solvent Accessible Surface Area)](#7-cálculo-de-sasa-solvent-accessible-surface-area)
8. [Todas as Visualizações](#8-todas-as-visualizações)

---

## 1. Download e Limpeza do PDB

### 1.1 Download do PDB 1TUP (apenas cadeia A, sem heteroátomos)

```bash
python ../biohub.py fetchpdb 1TUP --chains A --protein-only --output 1TUP_clean.pdb
```

**O que este comando faz:**
- Baixa o PDB 1TUP do RCSB
- Filtra apenas a cadeia A
- Remove íons, água e heteroátomos
- Salva como `1TUP_clean.pdb`

**Saída esperada:**
```
Fetching PDB 1TUP from RCSB PDB...
PDB downloaded successfully.
Filtering PDB content...
Filtered PDB saved to: 1TUP_clean.pdb
```

[ADD IMAGEM AQUI]

### 1.2 Informações do PDB baixado

```bash
head -20 1TUP_clean.pdb
```

[ADD IMAGEM AQUI]

---

## 2. Conversão PDB para FASTA

### 2.1 Converter o PDB limpo para FASTA

```bash
python ../biohub.py fasta 1TUP_clean.pdb --output 1TUP.fasta
```

**O que este comando faz:**
- Extrai a sequência de aminoácidos do arquivo PDB
- Converte códigos de 3 letras (ALA, GLY, etc.) para 1 letra (A, G, etc.)
- Salva no formato FASTA

**Saída esperada:**
```
FASTA file written to: 1TUP.fasta
```

### 2.2 Visualizar o arquivo FASTA gerado

```bash
cat 1TUP.fasta
```

**Exemplo de saída:**
```
>1TUP Chain A
MTAMEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDLLLPQDVEEFFEGPSEALRVSGAPAAQDPVTETPGPVAPAPATPWPLSSFVPSQKTYQGNYGFHLGFLQSGTAKSVMCTYSPPLNKLFI...
```

[ADD IMAGEM AQUI]

---

## 3. Análises com FASTA

### 3.1 Propriedades Físico-Químicas

```bash
python ../biohub.py physchem MTAMEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDLLLPQDVEEFFEGPSEALRVSGAPAAQDPVTETPGPVAPAPATPWPLSSFVPSQKTYQGNYGFHLGFLQSGTAKSVMCTYSPPLNKLFI --output 1TUP_properties.csv --plot-composition 1TUP_composition.png --plot-treemap 1TUP_treemap.png --plot-hydro 1TUP_hydropathy.png --window 9
```

**NOTA:** Cole a sequência completa extraída do arquivo FASTA no lugar da sequência de exemplo acima.

**O que este comando calcula:**
- Composição de aminoácidos
- Peso molecular
- Ponto isoelétrico (pI)
- Carga líquida em pH 7.0
- Coeficiente de extinção molar
- Índice de instabilidade
- Índice alifático
- GRAVY (Grand Average of Hydropathy)

**Arquivos gerados:**
- `1TUP_properties.csv` - Dados em formato CSV
- `1TUP_composition.png` - Gráfico de composição de aminoácidos
- `1TUP_treemap.png` - Treemap da composição
- `1TUP_hydropathy.png` - Perfil de hidropaticidade

#### Saída do CSV

```bash
cat 1TUP_properties.csv
```

[ADD IMAGEM AQUI]

#### Gráfico de Composição de Aminoácidos

![Composição de Aminoácidos](1TUP_composition.png)

[ADD IMAGEM AQUI]

#### Treemap da Composição

![Treemap](1TUP_treemap.png)

[ADD IMAGEM AQUI]

#### Perfil de Hidropaticidade

![Hidropaticidade](1TUP_hydropathy.png)

[ADD IMAGEM AQUI]

---

## 4. Conversão CSV para FASTA

Esta seção demonstra a conversão de dados CSV para formato FASTA.

### 4.1 Preparar um arquivo CSV de exemplo

**Crie um arquivo chamado `sequences.csv`:**

```csv
id,sequence,description
seq1,MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM,KRAS Protein
seq2,MTAMEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDLLLPQDVEEFFEGPSEALRVSGAPAAQDPVTETPGPVAPAPATPWPLSSFVPSQKTYQGNYGFHLGFLQSGTAKSVMCTYSPPLNKLFI,p53 Fragment
```

[ADD IMAGEM AQUI - CONTEÚDO DO CSV]

### 4.2 Converter CSV para FASTA

```bash
python ../biohub.py csv2fasta sequences.csv --output sequences.fasta --id-col id --seq-col sequence --header
```

**O que este comando faz:**
- Lê o arquivo CSV
- Extrai as colunas de ID e sequência
- Converte para formato FASTA

**Saída esperada:**
```
FASTA file created from CSV: sequences.fasta
```

### 4.3 Visualizar o FASTA gerado

```bash
cat sequences.fasta
```

**Exemplo de saída:**
```
>seq1 KRAS Protein
MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM
>seq2 p53 Fragment
MTAMEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDLLLPQDVEEFFEGPSEALRVSGAPAAQDPVTETPGPVAPAPATPWPLSSFVPSQKTYQGNYGFHLGFLQSGTAKSVMCTYSPPLNKLFI
```

[ADD IMAGEM AQUI]

### 4.4 Análises com o FASTA convertido

Agora você pode executar **todas as análises da seção 3** com as sequências convertidas do CSV:

**Para seq1 (KRAS):**
```bash
python ../biohub.py physchem MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM --output seq1_properties.csv --plot-composition seq1_composition.png --plot-treemap seq1_treemap.png --plot-hydro seq1_hydropathy.png --window 9
```

**Para seq2 (p53):**
```bash
python ../biohub.py physchem MTAMEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDLLLPQDVEEFFEGPSEALRVSGAPAAQDPVTETPGPVAPAPATPWPLSSFVPSQKTYQGNYGFHLGFLQSGTAKSVMCTYSPPLNKLFI --output seq2_properties.csv --plot-composition seq2_composition.png --plot-treemap seq2_treemap.png --plot-hydro seq2_hydropathy.png --window 9
```

[ADD IMAGEM AQUI - RESULTADOS DAS ANÁLISES]

---

## 5. Análise de Contatos Intramoleculares

### 5.1 Calcular contatos com distância de 8.0 Å

```bash
python ../biohub.py contacts 1TUP_clean.pdb --threshold 8.0 --output 1TUP_contacts.csv --plot 1TUP_contact_map.png
```

**O que este comando faz:**
- Identifica todos os pares de resíduos com distância entre Cα ≤ 8.0 Å
- Gera matriz de contatos
- Cria visualizações

**Arquivos gerados:**
- `1TUP_contacts.csv` - Lista de contatos com distâncias
- `1TUP_contact_map.png` - Mapa de contatos

#### Saída do CSV de Contatos

```bash
head -20 1TUP_contacts.csv
```

**Exemplo:**
```
Residue1,Residue2,Distance
M1,T2,3.82
M1,A3,5.45
T2,A3,3.81
T2,M4,3.79
...
```

[ADD IMAGEM AQUI]

#### Mapa de Contatos

![Mapa de Contatos](1TUP_contact_map.png)

[ADD IMAGEM AQUI]

---

## 6. Análise de Hidrofobicidade

### 6.1 Predizer a exposição ao solvente baseada em hidrofobicidade

```bash
python ../biohub.py hydrophoby 1TUP_clean.pdb --output 1TUP_hydrophoby.csv --write-pdb 1TUP_hydrophoby.pdb --pymol 1TUP_hydrophoby.pse --plot-hydrophoby 1TUP_hydrophoby_profile.png
```

**O que este comando faz:**
- Prediz exposição ao solvente baseado em hidrofobicidade (Kyte-Doolittle)
- Calcula hidrofobicidade por átomo
- Classifica como: hidrofóbico ou hidrofílico

**Arquivos gerados:**
- `1TUP_hydrophoby.csv` - Dados de hidrofobicidade por resíduo
- `1TUP_hydrophoby.pdb` - PDB com B-factors ajustados
- `1TUP_hydrophoby.pse` - Sessão PyMOL colorida por hidrofobicidade
- `1TUP_hydrophoby_profile.png` - Gráfico de perfil de hidrofobicidade

#### Saída do CSV de Hidrofobicidade

```bash
head -20 1TUP_hydrophoby.csv
```

**Exemplo:**
```
Chain,ResNum,ResName,AtomNum,AtomName,Hydrophobicity
A,1,MET,1,N,1.9
A,2,THR,2,CA,-0.7
A,3,ALA,3,C,1.8
A,4,MET,4,CB,1.9
...
```

[ADD IMAGEM AQUI]

#### Perfil de Hidrofobicidade

![Perfil Hidrofobicidade](1TUP_hydrophoby_profile.png)

[ADD IMAGEM AQUI]

#### Visualização PyMOL da Hidrofobicidade

```bash
pymol 1TUP_hydrophoby.pse
```

**Esquema de cores:**
- **Azul**: Resíduos hidrofílicos (valores negativos)
- **Branco**: Neutros
- **Vermelho**: Resíduos hidrofóbicos (valores positivos)

[ADD IMAGEM AQUI - SCREENSHOT DO PYMOL]

---

## 7. Cálculo de SASA (Solvent Accessible Surface Area)

### 7.1 Calcular SASA com 500 pontos

```bash
python ../biohub.py sasa 1TUP_clean.pdb --num-points 500 --output 1TUP_sasa.csv --write-pdb 1TUP_sasa.pdb --pymol 1TUP_sasa.pse --plot-profile 1TUP_sasa_profile.png
```

**O que este comando faz:**
- Calcula SASA usando algoritmo de Shrake-Rupley
- Usa 500 pontos distribuídos em uma esfera ao redor de cada átomo
- Gera perfil de SASA por resíduo

**Arquivos gerados:**
- `1TUP_sasa.csv` - SASA por resíduo em Ų
- `1TUP_sasa.pdb` - PDB com B-factors = SASA
- `1TUP_sasa.pse` - Sessão PyMOL colorida por SASA
- `1TUP_sasa_profile.png` - Gráfico do perfil de SASA

#### Saída do CSV de SASA

```bash
head -20 1TUP_sasa.csv
```

**Exemplo:**
```
Residue,SASA
M1,145.32
T2,98.45
A3,12.67
M4,35.89
E5,125.43
...
```

[ADD IMAGEM AQUI]

#### Perfil de SASA ao longo da sequência

![Perfil SASA](1TUP_sasa_profile.png)

[ADD IMAGEM AQUI]

#### Visualização PyMOL do SASA

```bash
pymol 1TUP_sasa.pse
```

**Esquema de cores (gradiente):**
- **Azul**: SASA baixo (resíduos enterrados)
- **Verde**: SASA intermediário
- **Vermelho**: SASA alto (resíduos expostos)

[ADD IMAGEM AQUI - SCREENSHOT DO PYMOL]

---

## 8. Todas as Visualizações

### 8.1 Resumo de todos os arquivos gerados

```bash
ls -lh
```

**Arquivos esperados:**

#### Arquivos de Dados
- `1TUP_clean.pdb` - PDB limpo (apenas cadeia A, sem heteroátomos)
- `1TUP.fasta` - Sequência em formato FASTA
- `sequences.csv` - CSV de exemplo com sequências
- `sequences.fasta` - FASTA convertido do CSV

#### CSVs de Análises
- `1TUP_properties.csv` - Propriedades físico-químicas
- `1TUP_contacts.csv` - Contatos intramoleculares
- `1TUP_hydrophoby.csv` - Predição de hidrofobicidade
- `1TUP_sasa.csv` - SASA por resíduo
- `seq1_properties.csv` - Propriedades da seq1 (KRAS)
- `seq2_properties.csv` - Propriedades da seq2 (p53)

#### Gráficos PNG
- `1TUP_composition.png` - Composição de aminoácidos
- `1TUP_treemap.png` - Treemap da composição
- `1TUP_hydropathy.png` - Perfil de hidropaticidade
- `1TUP_contact_map.png` - Mapa de contatos
- `1TUP_hydrophoby_profile.png` - Perfil de hidrofobicidade
- `1TUP_sasa_profile.png` - Perfil de SASA
- `seq1_composition.png`, `seq1_treemap.png`, `seq1_hydropathy.png`
- `seq2_composition.png`, `seq2_treemap.png`, `seq2_hydropathy.png`

#### Sessões PyMOL (.pse)
- `1TUP_hydrophoby.pse` - Visualização de hidrofobicidade
- `1TUP_sasa.pse` - Visualização de SASA

#### PDBs Modificados
- `1TUP_hydrophoby.pdb` - PDB com B-factors de hidrofobicidade
- `1TUP_sasa.pdb` - PDB com B-factors de SASA

[ADD IMAGEM AQUI - LISTAGEM COMPLETA]

---

## Checklist de Demonstração Completa

- [ ] Download do PDB 1TUP limpo (cadeia A, sem heteroátomos)
- [ ] Conversão PDB para FASTA
- [ ] Análise de propriedades físico-químicas do 1TUP (CSV + 3 gráficos)
- [ ] Conversão CSV para FASTA
- [ ] Análise de propriedades físico-químicas das sequências do CSV
- [ ] Análise de contatos intramoleculares (CSV + mapa)
- [ ] Predição de hidrofobicidade (CSV + PDB + PyMOL)
- [ ] Cálculo de SASA (CSV + gráfico + PDB + PyMOL)

---

## Notas Importantes

1. **Todos os comandos** devem ser executados dentro da pasta `demoday-biohub/`
2. O script principal está em `../biohub.py` (um nível acima)
3. Todas as análises estruturais usam o `1TUP_clean.pdb` (PDB limpo)
4. Para comandos `physchem`, cole a sequência completa extraída do FASTA
5. **PyMOL** deve estar instalado para abrir os arquivos `.pse`

---

## Comandos Sequenciais (Cópia Rápida)

```bash
# 1. Download e limpeza
python ../biohub.py fetchpdb 1TUP --chains A --protein-only --output 1TUP_clean.pdb

# 2. PDB para FASTA
python ../biohub.py fasta 1TUP_clean.pdb --output 1TUP.fasta

# 3. Propriedades físico-químicas (cole a sequência completa do 1TUP.fasta)
python ../biohub.py physchem SEQUENCIA_COMPLETA_AQUI --output 1TUP_properties.csv --plot-composition 1TUP_composition.png --plot-treemap 1TUP_treemap.png --plot-hydro 1TUP_hydropathy.png --window 9

# 4. Criar CSV de exemplo (cole o conteúdo manualmente em sequences.csv)

# 5. CSV para FASTA
python ../biohub.py csv2fasta sequences.csv --output sequences.fasta --id-col id --seq-col sequence --header

# 6. Análises das sequências do CSV
python ../biohub.py physchem SEQUENCIA_SEQ1_AQUI --output seq1_properties.csv --plot-composition seq1_composition.png --plot-treemap seq1_treemap.png --plot-hydro seq1_hydropathy.png --window 9
python ../biohub.py physchem SEQUENCIA_SEQ2_AQUI --output seq2_properties.csv --plot-composition seq2_composition.png --plot-treemap seq2_treemap.png --plot-hydro seq2_hydropathy.png --window 9

# 7. Contatos intramoleculares
python ../biohub.py contacts 1TUP_clean.pdb --threshold 8.0 --output 1TUP_contacts.csv --plot 1TUP_contact_map.png

# 8. Hidrofobicidade
python ../biohub.py hydrophoby 1TUP_clean.pdb --output 1TUP_hydrophoby.csv --write-pdb 1TUP_hydrophoby.pdb --pymol 1TUP_hydrophoby.pse --plot-hydrophoby 1TUP_hydrophoby_profile.png

# 9. SASA
python ../biohub.py sasa 1TUP_clean.pdb --num-points 500 --output 1TUP_sasa.csv --write-pdb 1TUP_sasa.pdb --pymol 1TUP_sasa.pse --plot-profile 1TUP_sasa_profile.png
```

---

## Contato

Para dúvidas ou sugestões sobre o BioHub, entre em contato.

**BioHub** - Uma plataforma completa para análise de sequências e estruturas de proteínas
