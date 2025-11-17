# RELATORIO DE ANALISE - ESTRUTURA 9VGK

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
python3 biohub.py fetchpdb 9VGK --chains A --protein-only -o 9VGK_clean.pdb
```

**O que este comando faz:**
- Baixa o PDB 9VGK do RCSB
- Seleciona apenas a cadeia A (--chains A)
- Remove agua, ligantes e heteroatomos (--protein-only)
- Salva em arquivo limpo: 9VGK_clean.pdb

**Resultado:**
- 2998 atomos mantidos
- 9498 atomos removidos (agua + ligantes + outras cadeias)

---

### 2. Extrair Sequencia FASTA
```bash
python3 biohub.py fasta 9VGK_clean.pdb -o 9VGK.fasta
```

**Sequencia extraida (388 residuos):**
```
>sequence_from_9VGK_clean.pdb
YNVDENGYYGEFGGAYIPEILHKCVEDLQNNYLKILESPDFQKEYDQLLRDYVGRPSPLYL
AKRLSEKYGCKIYLKREDLNHTGAHAINNTIGQILLARRMGKTRIIAETGAGQHGVATATA
VCALMNMECIVYMGKTDVERQHVNVQKMEMLGATVVPVTSGNMTLKDATNEAIRDWCCHPS
DTYYIIGSTVGPHPYPDMVARLQSVISKEIKKQLQEKEGRDYPDYLIACVGGGSNAAGTI
YHYIDDERVKIVLAEAGGKGIDSGMTAATIHLGKMGIIHGSKTLVMQNEDGQIEEPYSISA
GLDYPGIGPMHANLAKQKRAQVLAINDDEALNAAFELTRLEGIIPALESAHALAALEKVK
FKPEDVVVLTLSGRGDKDMETYLKY
```

---

### 3. Analise Fisico-Quimica (physchem)
```bash
python3 biohub.py physchem "SEQUENCIA_AQUI" -o 9VGK_physchem.csv \
  --plot-treemap 9VGK_treemap.png \
  --plot-composition 9VGK_composition.png \
  --plot-hydro 9VGK_hydro_profile.png
```

**Resultados:**

| Propriedade | Valor |
|-------------|-------|
| Comprimento | 388 residuos |
| Peso Molecular | 42,776.62 Da (~42.8 kDa) |
| Ponto Isoeletrico (pI) | 5.93 |
| GRAVY (Hidropaticidade) | -0.282 (hidrofilico) |
| Indice Alifatico | 90.03 |
| Indice de Instabilidade | 36.09 (Estavel) |
| Meia-Vida (E. coli) | 2-30 min |
| Residuos Acidos (Asp+Glu) | 51 |
| Residuos Basicos (Arg+Lys+His) | 53 |
| Residuos Polares | 94 |
| Residuos Apolares | 190 |

**Interpretacao:**
- Proteina media (~43 kDa)
- Levemente acida (pI 5.93, balanco entre acidos e basicos)
- Hidrofilica (GRAVY negativo)
- Estavel (indice de instabilidade < 40)
- Meia-vida curta devido ao residuo N-terminal (Tyr)

---

### 4. Analise de Hidrofobicidade (hydrophoby)
```bash
python3 biohub.py hydrophoby 9VGK_clean.pdb -o 9VGK_hydrophoby.csv \
  --plot-hydrophoby 9VGK_hydrophoby_plot.png
```

**Resultado:**
- 2998 atomos analisados
- Dados salvos em CSV com hidrofobicidade por atomo (escala Kyte-Doolittle)
- Grafico de perfil de hidrofobicidade gerado

---

### 5. Analise de Contatos Intramoleculares (contacts)
```bash
python3 biohub.py contacts 9VGK_clean.pdb -o 9VGK_contacts.csv \
  --plot 9VGK_contact_map.png
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
python3 biohub.py sasa 9VGK_clean.pdb --num-points 100 \
  -o 9VGK_sasa.csv --plot-profile 9VGK_sasa_plot.png
```

**Resultado:**
- SASA Total da Molecula: 15210.45 A²
- 2998 atomos analisados
- Tempo de execucao: ~7 minutos (com 100 pontos)
- Dados salvos em CSV
- Grafico de perfil SASA gerado

**Nota sobre Performance:**
- Com 960 pontos (padrão): execucao extremamente lenta (>20 minutos ou trava)
- Com 100 pontos: execucao aceitavel (~7 minutos)
- Complexidade: O(N² × P) onde N = atomos e P = pontos na esfera
- Para 2998 atomos com 960 pontos: ~8.6 bilhoes de calculos
- Para 2998 atomos com 100 pontos: ~900 milhoes de calculos
- Proteina quase 2x maior que 3E9C, tempo de execucao ~5x maior

---

## INFORMACOES DA ESTRUTURA

**PDB ID:** 9VGK
**Titulo:** Ancestral L-tryptophan synthase beta-subunit 1 complex with TRP-PLP
**Data de Deposito:** 14-JUN-25
**Classificacao:** BIOSYNTHETIC PROTEIN
**Organismo:** Synthetic construct (proteina ancestral reconstruida)

---

## ARQUIVOS GERADOS

### Dados (CSV)
```
issue_report/9VGK/
├── 9VGK_clean.pdb              # PDB limpo (cadeia A, apenas proteina) - 305KB
├── 9VGK.fasta                  # Sequencia em formato FASTA
├── 9VGK_physchem.csv           # Propriedades fisico-quimicas
├── 9VGK_hydrophoby.csv         # Hidrofobicidade por atomo - 73KB
├── 9VGK_contacts.csv           # Contatos intramoleculares - 71KB
└── 9VGK_sasa.csv               # SASA por atomo - 68KB
```

### Graficos (PNG)
```
├── 9VGK_treemap.png            # Treemap de composicao de aminoacidos - 244KB
├── 9VGK_composition.png        # Grafico de barras de composicao - 163KB
├── 9VGK_hydro_profile.png      # Perfil de hidrofobicidade Kyte-Doolittle - 419KB
├── 9VGK_hydrophoby_plot.png    # Perfil de hidrofobicidade por residuo - 947KB
├── 9VGK_contact_map.png        # Mapa de contatos intramoleculares - 328KB
└── 9VGK_sasa_plot.png          # Perfil de SASA por residuo - 1.0MB
```

### Relatorio
```
└── REPORT_9VGK.md              # Este relatorio completo
```

---

## RESUMO PARA O USUARIO

1. **Use sempre `fetchpdb` com `--protein-only` e `--chains`** para limpar o PDB antes de analises
2. **O erro KeyError foi causado** por caracteres invalidos na sequencia extraida de PDB nao tratado
3. **Comando correto:**
   ```bash
   python3 biohub.py fetchpdb 9VGK --chains A --protein-only -o 9VGK_clean.pdb
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
