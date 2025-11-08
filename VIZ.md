# Exemplos de Uso das Visualizações do BioHub

Este documento contém exemplos práticos de como gerar visualizações gráficas com o BioHub.

## Instalação das Dependências

Antes de começar, instale as dependências de visualização:

```bash
pip install -r requirements.txt
```

Ou manualmente:

```bash
pip install matplotlib numpy squarify
```

## Exemplos por Tipo de Análise

### 1. Análise de Composição de Aminoácidos

#### Treemap (Visualização Moderna)

```bash
python3 biohub.py physchem "MKTAYIAKQRQISFVKSHFSRQ" --plot-treemap composicao_treemap.png
```

**Ideal para:**
- Comparação visual rápida de composição
- Identificação de aminoácidos predominantes
- Apresentações e publicações

#### Gráfico de Barras (Visualização Clássica)

```bash
python3 biohub.py physchem "MKTAYIAKQRQISFVKSHFSRQ" --plot-composition composicao_barras.png
```

**Ideal para:**
- Análise quantitativa detalhada
- Comparação com outros trabalhos
- Publicações científicas tradicionais

### 2. Análise de Hidrofobicidade

#### Perfil Kyte-Doolittle Padrão

```bash
python3 biohub.py physchem "MKTAYIAKQRQISFVKSHFSRQ..." --plot-hydro hydro.png
```

Usa janela de 9 resíduos (padrão Kyte-Doolittle original).

#### Para Proteínas de Membrana

```bash
python3 biohub.py physchem "SEQUENCE..." --plot-hydro hydro_membrane.png --window 19
```

Janelas maiores (19-21) são melhores para identificar domínios transmembrana.

**Interpretação:**
- Picos > 1.6: possíveis hélices transmembrana
- Regiões hidrofílicas: loops extracelulares/citoplasmáticos

### 3. Análise de SASA (Acessibilidade ao Solvente)

#### Perfil de SASA

```bash
python3 biohub.py sasa proteina.pdb --plot-profile sasa_profile.png
```

**Com alta precisão:**

```bash
python3 biohub.py sasa proteina.pdb --plot-profile sasa_profile.png --num-points 2000
```

Mais pontos = maior precisão (mas mais lento).

**Interpretação:**
- Picos altos (>20 Ų): resíduos expostos (potenciais epitopos)
- Vales baixos (<5 Ų): resíduos enterrados (núcleo hidrofóbico)

### 4. Análise de Contatos Intramoleculares

#### Mapa de Contatos Padrão

```bash
python3 biohub.py contacts proteina.pdb --plot contact_map.png
```

#### Com Threshold Customizado

```bash
python3 biohub.py contacts proteina.pdb --plot contact_map.png -t 10.0
```

Aumentar threshold (ex: 10.0 Å) mostra mais contatos de longo alcance.

**Interpretação:**
- Padrões diagonais: hélices alfa
- Padrões paralelos: folhas beta
- Blocos fora da diagonal: domínios estruturais

## Análises Combinadas

### Análise Completa de Sequência

```bash
# Gera todos os gráficos de uma vez + dados CSV
python3 biohub.py physchem "MKTAYIAKQRQISFVK..." \
  --plot-treemap treemap.png \
  --plot-composition barras.png \
  --plot-hydro hydro.png \
  --window 11 \
  -o propriedades.csv
```

### Análise Estrutural Completa

```bash
# Combina SASA com PyMOL e gráficos
python3 biohub.py sasa proteina.pdb \
  --plot-profile sasa_profile.png \
  --write-pdb proteina_sasa.pdb \
  --pymol visualizacao.pse \
  -o sasa_data.csv
```

```bash
# Mapa de contatos + CSV
python3 biohub.py contacts proteina.pdb \
  --plot contact_map.png \
  -t 8.0 \
  -o contatos.csv
```

## Casos de Uso Específicos

### Caso 1: Comparação de Proteínas Homólogas

```bash
# Proteína 1
python3 biohub.py physchem "$(python3 biohub.py fasta prot1.pdb)" \
  --plot-treemap prot1_comp.png

# Proteína 2
python3 biohub.py physchem "$(python3 biohub.py fasta prot2.pdb)" \
  --plot-treemap prot2_comp.png
```

Compare visualmente os treemaps para identificar diferenças na composição.

### Caso 2: Predição de Epitopos

```bash
# Identifique regiões expostas
python3 biohub.py sasa anticorpo.pdb --plot-profile epitopos.png

# Combine com hidrofobicidade
python3 biohub.py physchem "SEQUENCE" --plot-hydro hydro_epitopos.png
```

Epitopos geralmente são:
- Alta exposição (SASA alto)
- Hidrofílicos (hidrofobicidade negativa)

### Caso 3: Análise de Proteína de Membrana

```bash
# Sequência completa
SEQ="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK..."

# Perfil com janela grande
python3 biohub.py physchem "$SEQ" \
  --plot-hydro membrane_topology.png \
  --window 21

# Composição (proteínas de membrana são ricas em hidrofóbicos)
python3 biohub.py physchem "$SEQ" --plot-treemap membrane_comp.png
```

### Caso 4: Otimização para Cristalização

```bash
# Identifique regiões desordenadas/flexíveis
python3 biohub.py sasa proteina.pdb \
  --plot-profile flexibility.png \
  --num-points 1500

# Mapa de contatos para ver domínios estruturados
python3 biohub.py contacts proteina.pdb \
  --plot domains.png \
  -t 8.0
```

Regiões com SASA muito variável podem ser flexíveis → candidatas a truncagem.

## Formatos de Saída

### PNG (Padrão - Recomendado)

```bash
--plot-treemap figura.png
```

- Tamanho razoável
- Boa qualidade (300 DPI)
- Suportado por todos os programas

### PDF (Vetorial - Para Publicações)

```bash
--plot-treemap figura.pdf
```

- Escalável sem perda de qualidade
- Ideal para revistas científicas
- Pode ser editado em Illustrator/Inkscape

### SVG (Vetorial - Editável)

```bash
--plot-treemap figura.svg
```

- Totalmente editável
- Bom para apresentações
- Pode abrir em navegadores

## Troubleshooting

### Erro: "Módulo de visualização não disponível"

```bash
pip install matplotlib numpy squarify
```

### Erro: "squarify não instalado"

Treemap usará visualização alternativa. Para treemap completo:

```bash
pip install squarify
```

### Verificar Dependências

```bash
python3 -c "import matplotlib, numpy, squarify; print('OK!')"
```

### Gráficos Muito Grandes

Se os arquivos PNG ficarem muito grandes, ajuste a resolução:

```python
# Edite biohub_viz.py e mude:
plt.savefig(output_file, dpi=150)  # Ao invés de 300
```

### Fontes Muito Pequenas

Para apresentações, aumente as fontes editando `biohub_viz.py`:

```python
plt.rcParams['font.size'] = 14  # Adicione no início das funções
```

## Dicas de Qualidade

1. **Para publicações:** Use PDF ou SVG
2. **Para apresentações:** Use PNG com DPI 300
3. **Para web:** Use PNG com DPI 150
4. **Combine múltiplos gráficos:** Use Inkscape ou PowerPoint para layout

## Exemplos com Dados Reais

### Hemoglobina (4HHB)

```bash
# Baixa estrutura
python3 biohub.py fetchpdb 4HHB --protein-only -o hemoglobina.pdb

# Análise de contatos
python3 biohub.py contacts hemoglobina.pdb --plot hb_contacts.png

# SASA
python3 biohub.py sasa hemoglobina.pdb --plot-profile hb_sasa.png --num-points 500
```

### Crambin (1CRN) - Proteína Pequena

```bash
# Baixa estrutura
python3 biohub.py fetchpdb 1CRN --protein-only -o crambin.pdb

# Extrai sequência
SEQ=$(python3 biohub.py fasta crambin.pdb)

# Análise completa
python3 biohub.py physchem "$SEQ" \
  --plot-treemap crambin_comp.png \
  --plot-hydro crambin_hydro.png \
  --window 7
```

## Referências

- Kyte & Doolittle (1982) - Perfil de hidrofobicidade
- Lee & Richards (1971) - Método SASA
- Shrake & Rupley (1973) - Algoritmo SASA

## Suporte

Para reportar problemas ou sugerir melhorias nas visualizações, abra um issue no repositório do projeto.
