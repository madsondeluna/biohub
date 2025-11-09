# BioHub - Exemplos de Visualizações

Esta pasta contém exemplos de todas as visualizações geradas pelo BioHub usando a estrutura da hemoglobina (PDB: 4HHB).

## Arquivos Gerados

### Arquivos de Sequência
- **4HHB.fasta** - Sequência FASTA extraída da cadeia A do PDB 4HHB

### Análises Físico-Químicas (physchem)

#### Dados
- **4HHB_physchem.csv** - Propriedades físico-químicas calculadas (peso molecular, pI, GRAVY, índice alifático, etc.)

#### Visualizações
- **4HHB_treemap.png** - Treemap hierárquico da composição de aminoácidos agrupados por propriedade química
- **4HHB_composition.png** - Gráfico de barras da composição de aminoácidos ordenado por frequência
- **4HHB_hydropathy.png** - Perfil de hidrofobicidade Kyte-Doolittle com janela deslizante (window=9)

### Análise de Contatos Intramoleculares (contacts)

#### Dados
- **4HHB_contacts.csv** - Lista de todos os contatos intramoleculares (Cα-Cα) com distância ≤ 8.0 Å

#### Visualizações
- **4HHB_contact_map.png** - Mapa de contatos (heatmap) mostrando interações entre resíduos

### Análise de SASA (sasa)

#### Dados
- **4HHB_sasa.csv** - SASA calculado por átomo usando algoritmo de Shrake-Rupley (300 pontos)
- **4HHB_sasa.pdb** - Estrutura PDB anotada com SASA médio por resíduo no B-factor

#### Visualizações
- **4HHB_sasa_profile.png** - Perfil de SASA por resíduo ao longo da sequência
- **4HHB_sasa.pml** - Script PyMOL para visualização 3D da estrutura com SASA colorido

#### Visualização PyMOL
Para abrir a visualização 3D no PyMOL:
```bash
pymol 4HHB_sasa.pml
```

## Informações da Estrutura 4HHB

- **Proteína**: Hemoglobina humana
- **Organismo**: *Homo sapiens*
- **Resolução**: 1.74 Å
- **Método**: Difração de raios-X
- **Cadeias**: A, B (alfa e beta globinas)
- **Resíduos analisados**: 141 (cadeia A)
- **Total de átomos**: 1069

## Estatísticas das Análises

### SASA
- **SASA Total**: 7911.28 Ų
- **Resíduos enterrados** (SASA < 0.5 Ų): 9
- **Resíduos expostos** (SASA ≥ 0.5 Ų): 132
- **Range de SASA**: 0.00 - 23.09 Ų

### Contatos
- **Threshold**: 8.0 Å
- **Total de contatos identificados**: Veja o arquivo CSV para detalhes

## Qualidade das Imagens

Todas as visualizações foram geradas com:
- **Resolução**: 300 DPI (qualidade para publicação)
- **Formato**: PNG com fundo branco
- **Paleta de cores**: Científica e amigável para daltônicos

## Comandos Utilizados

```bash
# Extrair sequência FASTA
python3 biohub.py fasta examples-pymol/4HHB.pdb -o viz-exemplos/4HHB.fasta

# Análises físico-químicas com visualizações
python3 biohub.py physchem "SEQUENCIA..." \
  --plot-treemap viz-exemplos/4HHB_treemap.png \
  --plot-composition viz-exemplos/4HHB_composition.png \
  --plot-hydro viz-exemplos/4HHB_hydropathy.png \
  --window 9 \
  -o viz-exemplos/4HHB_physchem.csv

# Mapa de contatos
python3 biohub.py contacts examples-pymol/4HHB.pdb \
  --plot viz-exemplos/4HHB_contact_map.png \
  -t 8.0 \
  -o viz-exemplos/4HHB_contacts.csv

# SASA com perfil e visualização PyMOL
python3 biohub.py sasa examples-pymol/4HHB.pdb \
  --plot-profile viz-exemplos/4HHB_sasa_profile.png \
  --write-pdb viz-exemplos/4HHB_sasa.pdb \
  --pymol viz-exemplos/4HHB_sasa.pse \
  -o viz-exemplos/4HHB_sasa.csv \
  --num-points 300
```

## Uso Educacional

Estes exemplos demonstram as capacidades do BioHub para análise e visualização de propriedades proteicas. Use-os como referência para:
- Entender os diferentes tipos de análises disponíveis
- Aprender a interpretar os gráficos gerados
- Verificar a instalação correta das dependências de visualização
- Criar seus próprios pipelines de análise

---

**Nota**: Algumas visualizações (como treemap) requerem a instalação de bibliotecas opcionais. Veja o arquivo `requirements.txt` na raiz do projeto.
