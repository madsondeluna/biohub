# Changelog - BioHub

## [0.1.3] - 2025-11-17

### Alterado
- **SASA: Redução do padrão de num-points de 960 para 200**
  - Motivo: O valor de 960 pontos causava cálculos extremamente lentos (>10-20 minutos) ou travamentos em estruturas médias a grandes
  - Complexidade do algoritmo: O(N² × P) onde N = número de átomos e P = pontos na esfera
  - Impacto da mudança:
    - Para estruturas pequenas (~200 resíduos, ~1600 átomos): redução de >10min para ~3-4min
    - Para estruturas médias (~400 resíduos, ~3000 átomos): redução de >20min para ~14min
    - Redução de aproximadamente 5x no número de cálculos (de bilhões para centenas de milhões)
  - Precisão: 200 pontos mantém precisão aceitável para análises exploratórias
  - Flexibilidade: Usuários que necessitarem maior precisão podem aumentar manualmente com `--num-points`
  - Localização da mudança: `biohub.py` linha 1055
  - Comentários explicativos adicionados no código (linhas 1051-1054)

### Observações
- Esta mudança resolve o issue reportado sobre SASA aparentemente "travando"
- O problema não era um bug, mas sim uma questão de performance com o valor padrão muito alto
- Versão atualizada de 0.1.2 para 0.1.3

## [0.1.2] - 2025-XX-XX

### Inicial
- Primeira versão documentada
- Módulos: fetchpdb, fasta, csv2fasta, physchem, contacts, hydrophoby, sasa, apbs
- Suporte a visualização com matplotlib, numpy, squarify
