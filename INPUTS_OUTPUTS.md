# Tabela Completa de Funções do BioHub

| Função | Tipo | Parâmetro | Descrição | Obrigatório | Padrão/Valores |
|--------|------|-----------|-----------|:-----------:|----------------|
| **fetchpdb** | INPUT | `PDB_ID` | ID de 4 caracteres do PDB (ex: 1A2B) | ✓ | - |
| | OUTPUT | `-o, --output` | Nome do arquivo de saída | ✗ | `[PDB_ID].pdb` |
| | FLAG | `--chains` | Cadeias a manter (separadas por vírgula) | ✗ | Todas |
| | FLAG | `--protein-only` | Remove água, ligantes e heteroátomos | ✗ | `False` |
| **fasta** | INPUT | `ARQUIVO_PDB` | Caminho para arquivo PDB | ✓ | - |
| | OUTPUT | `-o, --output` | Salva em arquivo FASTA | ✗ | stdout |
| **csv2fasta** | INPUT | `ARQUIVO_CSV` | Caminho para arquivo CSV | ✓ | - |
| | OUTPUT | `-o, --output` | Salva em arquivo FASTA | ✗ | stdout |
| | FLAG | `--id-col` | Coluna do identificador (nome ou índice) | ✗ | `0` |
| | FLAG | `--seq-col` | Coluna da sequência (nome ou índice) | ✗ | `1` |
| | FLAG | `--header` | Primeira linha é cabeçalho | ✗ | `False` |
| | FLAG | `--delimiter` | Delimitador do CSV | ✗ | `,` |
| **physchem** | INPUT | `SEQUENCIA` | Sequência de aminoácidos (1 letra) | ✓ | - |
| | OUTPUT | `-o, --output` | Salva resultados em CSV | ✗ | stdout |
| | PLOT | `--plot-treemap` | Gera treemap de composição (PNG) | ✗ | - |
| | PLOT | `--plot-composition` | Gera gráfico de barras de composição (PNG) | ✗ | - |
| | PLOT | `--plot-hydro` | Gera perfil de hidrofobicidade (PNG) | ✗ | - |
| | FLAG | `--window` | Tamanho da janela para hidrofobicidade | ✗ | `9` |
| **contacts** | INPUT | `ARQUIVO_PDB` | Caminho para arquivo PDB | ✓ | - |
| | OUTPUT | `-o, --output` | Salva resultados em CSV | ✗ | stdout |
| | FLAG | `-t, --threshold` | Distância máxima para contato (Å) | ✗ | `8.0` |
| | PLOT | `--plot` | Gera mapa de contatos (PNG) | ✗ | - |
| **hydrophoby** | INPUT | `ARQUIVO_PDB` | Caminho para arquivo PDB | ✓ | - |
| | OUTPUT | `-o, --output` | Salva resultados em CSV | ✗ | stdout |
| | OUTPUT | `--write-pdb` | Gera PDB com hidrofobicidade no B-factor | ✗ | - |
| | OUTPUT | `--pymol` | Gera sessão PyMOL (.pse + .pml) | ✗ | - |
| | PLOT | `--plot-hydrophoby` | Gera perfil de hidrofobicidade (PNG) | ✗ | - |
| **sasa** | INPUT | `ARQUIVO_PDB` | Caminho para arquivo PDB | ✓ | - |
| | OUTPUT | `-o, --output` | Salva resultados por átomo em CSV | ✗ | stdout |
| | OUTPUT | `--write-pdb` | Gera PDB com SASA no B-factor | ✗ | - |
| | OUTPUT | `--pymol` | Gera sessão PyMOL (.pse + .pml) | ✗ | - |
| | FLAG | `--probe-radius` | Raio da sonda do solvente (Å) | ✗ | `1.4` (água) |
| | FLAG | `--num-points` | Pontos na superfície de cada átomo | ✗ | `200` |
| | PLOT | `--plot-profile` | Gera perfil de SASA por resíduo (PNG) | ✗ | - |
| **apbs** | INPUT | `ARQUIVO_PDB` | Caminho para arquivo PDB | ✓ | - |
| | OUTPUT | - | Energia em kJ/mol (stdout) | - | - |
| | FLAG | `--no-cleanup` | Mantém arquivos temporários (PQR, apbs.in) | ✗ | `False` |

---

## Legendas

- **INPUT**: Arquivo ou dado de entrada obrigatório
- **OUTPUT**: Arquivo de saída opcional (sem especificar = stdout)
- **FLAG**: Parâmetro opcional que modifica comportamento
- **PLOT**: Gera arquivo de visualização (requer matplotlib/numpy)
- **✓**: Obrigatório | **✗**: Opcional

---

## Formatos de Arquivo Suportados

### Inputs Aceitos
- **PDB** (`.pdb`) - Protein Data Bank format
- **CSV** (`.csv`) - Valores separados por vírgula
- **Sequência** - String de aminoácidos (código de 1 letra)

### Outputs Gerados
- **FASTA** (`.fasta`) - Formato de sequência
- **CSV** (`.csv`) - Dados tabulares
- **PDB** (`.pdb`) - PDB anotado com B-factor
- **PNG** (`.png`) - Gráficos de visualização (300 DPI)
- **PSE** (`.pse`) - Sessão PyMOL
- **PML** (`.pml`) - Script PyMOL

---

## Dependências Opcionais

| Biblioteca | Comandos que Usam | Instalação |
|------------|-------------------|------------|
| **matplotlib** | Todos os `--plot-*` | `pip install matplotlib` |
| **numpy** | contacts, sasa, hydrophoby (plots) | `pip install numpy` |
| **squarify** | physchem `--plot-treemap` | `pip install squarify` |
| **PyMOL** | `--pymol` (sasa, hydrophoby) | Sistema-específico |
| **PDB2PQR** | apbs | Sistema-específico |
| **APBS** | apbs | Sistema-específico |

---

## Notas Importantes

1. **Performance SASA**: `--num-points` padrão reduzido para 200 (v0.1.3) para balancear precisão e velocidade
2. **PyMOL Sessions**: Se PyMOL não estiver no PATH, apenas o arquivo `.pml` será gerado
3. **Filtros PDB**: `fetchpdb` aplica filtros (`--chains`, `--protein-only`) após o download
4. **Output padrão**: Sem flag `-o`, a maioria dos comandos imprime em stdout (exceto plots e sessões PyMOL)
5. **CSV Encoding**: Todos os CSVs são gerados com encoding UTF-8
6. **APBS**: Comando em BETA, pode ser descontinuado em versões futuras
