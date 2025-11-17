# RELATORIO DE ANALISE - ESTRUTURA 3E9C

## PROBLEMA IDENTIFICADO NO ISSUE

O usuario executou os comandos incorretamente:
1. Nao usou o modulo `fetchpdb` para baixar e limpar o PDB
2. Nao especificou chains nem removeu heteroatomos
3. A sequencia continha caracteres nao-aminoacidicos (numeros) que causaram o erro no physchem

### Erro Reportado
```
KeyError: '9'
```

**Causa**: A sequencia extraida do PDB original continha caracteres invalidos (numeros de residuos, heteroatomos, etc.) que nao sao aminoacidos validos.

---

## COMANDOS CORRETOS PARA EXECUTAR

### 1. Baixar e Limpar o PDB
```bash
python3 biohub.py fetchpdb 3E9C --chains A --protein-only -o 3E9C_clean.pdb
```

**O que este comando faz:**
- Baixa o PDB 3E9C do RCSB
- Seleciona apenas a cadeia A (--chains A)
- Remove agua, ligantes e heteroatomos (--protein-only)
- Salva em arquivo limpo: 3E9C_clean.pdb

**Resultado:**
- 1597 atomos mantidos
- 2531 atomos removidos (agua + ligantes + outras cadeias)

---

### 2. Extrair Sequencia FASTA
```bash
python3 biohub.py fasta 3E9C_clean.pdb -o 3E9C.fasta
```

**Sequencia extraida (204 residuos):**
```
>sequence_from_3E9C_clean.pdb
MLTFALTIVRHGETDTPLSDTGHQQAAAAGRYLKDLHFTNVFVSNLQRAIQTAEIILGNNL
HSSATEMILDPLLRERGFGETLEQVKTRFKMFLKSLFQRMFEEHGSALSSADQPVIAGLAD
DGAQNVPVHALMVSHGAFIRISVRHLVEDLQCCLPAGLKMNQVFSPCPNTGISRFIFTIH
REESVLRATRIQGVFINRKDHL
```

---

### 3. Analise Fisico-Quimica (physchem)
```bash
python3 biohub.py physchem "SEQUENCIA_AQUI" -o 3E9C_physchem.csv \
  --plot-treemap 3E9C_treemap.png \
  --plot-composition 3E9C_composition.png \
  --plot-hydro 3E9C_hydro_profile.png
```

**Resultados:**

| Propriedade | Valor |
|-------------|-------|
| Comprimento | 204 residuos |
| Peso Molecular | 22,734.99 Da (~22.7 kDa) |
| Ponto Isoeletrico (pI) | 7.63 |
| GRAVY (Hidropaticidade) | -0.028 (neutro/levemente hidrofilico) |
| Indice Alifatico | 95.64 |
| Indice de Instabilidade | 29.86 (Estavel) |
| Meia-Vida (E. coli) | >10 horas |
| Residuos Acidos (Asp+Glu) | 20 |
| Residuos Basicos (Arg+Lys+His) | 30 |
| Residuos Polares | 49 |
| Residuos Apolares | 105 |

**Interpretacao:**
- Proteina pequena (~23 kDa)
- Levemente basica (pI 7.63, mais residuos basicos que acidos)
- Hidrofilica
- Estavel (indice de instabilidade < 40)
- Alta proporcao de residuos alifaticos

---

### 4. Analise de Hidrofobicidade (hydrophoby)
```bash
python3 biohub.py hydrophoby 3E9C_clean.pdb -o 3E9C_hydrophoby.csv \
  --plot-hydrophoby 3E9C_hydrophoby_plot.png
```

**Resultado:**
- 1597 atomos analisados
- Dados salvos em CSV com hidrofobicidade por atomo (escala Kyte-Doolittle)
- Grafico de perfil de hidrofobicidade gerado

---

### 5. Analise de Contatos Intramoleculares (contacts)
```bash
python3 biohub.py contacts 3E9C_clean.pdb -o 3E9C_contacts.csv \
  --plot 3E9C_contact_map.png
```

**Resultado:**
- Contatos intramoleculares calculados (threshold = 8.0 A)
- Dados salvos em CSV
- Mapa de contatos gerado

---

### 6. Analise SASA (Solvent Accessible Surface Area)

**SOLUCAO ENCONTRADA:**
O comando SASA com parametros padrão (960 pontos) é extremamente lento e parece travar.
**Solucao:** Reduzir o número de pontos para 100.

```bash
python3 biohub.py sasa 3E9C_clean.pdb --num-points 100 \
  -o 3E9C_sasa.csv --plot-profile 3E9C_sasa_plot.png
```

**Resultado:**
- SASA Total da Molecula: 9477.60 A²
- 1597 atomos analisados
- Tempo de execucao: ~1.5 minutos (com 100 pontos)
- Dados salvos em CSV
- Grafico de perfil SASA gerado

**Nota sobre Performance:**
- Com 960 pontos (padrão): execucao extremamente lenta (>10 minutos ou trava)
- Com 100 pontos: execucao aceitavel (~1.5 minutos)
- Complexidade: O(N² × P) onde N = atomos e P = pontos na esfera
- Para 1597 atomos com 960 pontos: ~2.4 bilhoes de calculos
- Para 1597 atomos com 100 pontos: ~255 milhoes de calculos

---

## INFORMACOES DA ESTRUTURA

**PDB ID:** 3E9C
**Titulo:** Structure of a tryptic core fragment of TIGAR from Danio rerio
**Data de Deposito:** 21-AUG-08
**Classificacao:** HYDROLASE
**Organismo:** Danio rerio (peixe-zebra)

---

## ARQUIVOS GERADOS

### Dados (CSV)
```
issue_report/3E9C/
├── 3E9C_clean.pdb              # PDB limpo (cadeia A, apenas proteina) - 180KB
├── 3E9C.fasta                  # Sequencia em formato FASTA
├── 3E9C_physchem.csv           # Propriedades fisico-quimicas
├── 3E9C_hydrophoby.csv         # Hidrofobicidade por atomo - 38KB
├── 3E9C_contacts.csv           # Contatos intramoleculares - 32KB
└── 3E9C_sasa.csv               # SASA por atomo - 36KB
```

### Graficos (PNG)
```
├── 3E9C_treemap.png            # Treemap de composicao de aminoacidos - 221KB
├── 3E9C_composition.png        # Grafico de barras de composicao - 150KB
├── 3E9C_hydro_profile.png      # Perfil de hidrofobicidade Kyte-Doolittle - 396KB
├── 3E9C_hydrophoby_plot.png    # Perfil de hidrofobicidade por residuo - 800KB
├── 3E9C_contact_map.png        # Mapa de contatos intramoleculares - 260KB
└── 3E9C_sasa_plot.png          # Perfil de SASA por residuo - 756KB
```

### Relatorio
```
└── REPORT_3E9C.md              # Este relatorio completo
```

---

## RESUMO PARA O USUARIO

1. **Use sempre `fetchpdb` com `--protein-only` e `--chains`** para limpar o PDB antes de analises
2. **O erro KeyError foi causado** por caracteres invalidos na sequencia extraida de PDB nao tratado
3. **Comando correto:**
   ```bash
   python3 biohub.py fetchpdb 3E9C --chains A --protein-only -o 3E9C_clean.pdb
   ```
4. **SASA apresenta bug confirmado** - executa indefinidamente sem completar
5. **Todos os outros modulos funcionam corretamente** quando o PDB esta limpo

---

## MELHORIAS NECESSARIAS NO BIOHUB

### 1. Tratamento de Erros Inadequado

**Problema Identificado:**
A aplicacao nao possui tratamento de erros robusto. Quando o usuario:
- Usa comandos incorretos
- Esquece flags obrigatorias
- Fornece dados invalidos (ex: sequencia com caracteres nao-aminoacidicos)

A aplicacao retorna apenas:
- Traceback completo do Python
- Numero da linha do codigo quebrado
- Mensagem de erro tecnica (ex: `KeyError: '9'`)

**Exemplo do Erro Atual:**
```
Traceback (most recent call last):
  File "/mnt/c/Users/leo/Documents/GitHub/biohub/biohub.py", line 1090, in <module>
    main()
  File "/mnt/c/Users/leo/Documents/GitHub/biohub/biohub.py", line 1079, in main
    command_functions[args.command](args)
  File "/mnt/c/Users/leo/Documents/GitHub/biohub/biohub.py", line 407, in calculate_physicochemical_properties
    instability_index = (10 / length) * sum(DIWV[sequence[i]][sequence[i+1]] for i in range(length - 1)) if length > 1 else 0
  File "/mnt/c/Users/leo/Documents/GitHub/biohub/biohub.py", line 407, in <genexpr>
    instability_index = (10 / length) * sum(DIWV[sequence[i]][sequence[i+1]] for i in range(length - 1)) if length > 1 else 0
KeyError: '9'
```

**Melhorias Recomendadas:**

1. **Validacao de Sequencia no physchem:**
   ```python
   # Adicionar antes do calculo:
   invalid_chars = set(sequence) - set(MOLECULAR_WEIGHT.keys())
   if invalid_chars:
       print(f"Erro: Sequencia contem caracteres invalidos: {invalid_chars}", file=sys.stderr)
       print(f"Apenas aminoacidos de letra unica (A-Z) sao aceitos.", file=sys.stderr)
       print(f"Dica: Use 'fetchpdb --protein-only --chains A' para limpar o PDB antes de extrair a sequencia.", file=sys.stderr)
       return
   ```

2. **Mensagens de Erro Amigaveis:**
   - Substituir tracebacks tecnicos por mensagens claras
   - Incluir sugestoes de correcao
   - Indicar qual parametro esta faltando ou incorreto

3. **Validacao de Entrada:**
   - Verificar se arquivos PDB existem antes de processar
   - Validar formato de sequencias antes de calculos
   - Checar se flags obrigatorias foram fornecidas

### 2. Bug Confirmado: SASA

**Status:** O modulo SASA executa indefinidamente sem completar o calculo, se o usuário não sinalizar a flag de output.

**Impacto:** Inviabiliza a analise de superficies acessiveis ao solvente, e o módulo não retorna resultados. 

**Prioridade:** Média-Alta. Deve ser corrigido para garantir funcionalidade completa, ou ao menos tratar o caso com uma mensagem de erro.

---

**Data do Relatorio:** 2025-11-17
**Ferramenta:** BioHub v0.1.2
