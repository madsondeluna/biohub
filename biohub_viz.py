#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioHub Visualization Module
Módulo de visualização para gráficos das análises bioinformáticas do BioHub.
Requer matplotlib, numpy e squarify como dependências opcionais.
"""

import sys

# Tentar importar dependências opcionais
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib
    matplotlib.use('Agg')  # Backend não-interativo (sem GUI)
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Aviso: matplotlib não disponível. Instale com: pip install matplotlib", file=sys.stderr)

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    print("Aviso: numpy não disponível. Instale com: pip install numpy", file=sys.stderr)

try:
    import squarify
    HAS_SQUARIFY = True
except ImportError:
    HAS_SQUARIFY = False


# ============================================================================
# CONSTANTES DE PADRONIZAÇÃO VISUAL
# ============================================================================
# Posicionamento padronizado para legendas e notas no lado direito
LEGEND_X_POSITION = 1.05   # Posição X padrão para todos os elementos à direita
LEGEND_SPACING = 0.03      # Espaçamento vertical entre blocos (legenda, nota, stats)

# Posição X especial para gráficos com colorbar (SASA)
LEGEND_X_WITH_COLORBAR = 1.18  # Deslocado para direita para evitar colorbar

# Posições Y calculadas com base no espaçamento
LEGEND_TOP_Y = 0.98        # Posição Y para legenda no topo
INFO_BOX_1_Y = 0.65        # Primeira caixa de info abaixo da legenda
INFO_BOX_2_Y = 0.40        # Segunda caixa de info (se houver)
INFO_BOX_3_Y = 0.15        # Terceira caixa de info (se houver)

# Estilo padronizado para legendas
LEGEND_FONT_SIZE = 10
LEGEND_STYLE = {
    'frameon': True,
    'shadow': True,
    'fancybox': True,
    'fontsize': LEGEND_FONT_SIZE
}

# Estilo padronizado para caixas de informação (stats e notas)
INFO_BOX_STYLE = {
    'boxstyle': 'round',
    'facecolor': 'wheat',
    'alpha': 0.9,
    'edgecolor': 'black',
    'linewidth': 1.0,
    'pad': 0.5
}

INFO_BOX_FONT_SIZE = 9
INFO_BOX_FAMILY = 'sans-serif'  # Mesma fonte do gráfico (não monospace)

# Margens para tight_layout (formato: [left, bottom, right, top])
LAYOUT_MARGINS_WITH_LEGEND = [0, 0, 0.85, 1]         # Para gráficos com legenda à direita
LAYOUT_MARGINS_WITH_COLORBAR = [0, 0, 0.78, 1]       # Para gráficos com colorbar
LAYOUT_MARGINS_SIMPLE = [0, 0, 1, 1]                 # Para gráficos sem legenda lateral


def check_dependencies(require_numpy=False, require_squarify=False):
    """
    Verifica se as dependências necessárias estão disponíveis.

    Args:
        require_numpy: Se True, verifica numpy
        require_squarify: Se True, verifica squarify

    Returns:
        True se todas as dependências estão disponíveis, False caso contrário
    """
    if not HAS_MATPLOTLIB:
        print("Erro: matplotlib não instalado. Instale com: pip install matplotlib", file=sys.stderr)
        return False

    if require_numpy and not HAS_NUMPY:
        print("Erro: numpy não instalado. Instale com: pip install numpy", file=sys.stderr)
        return False

    if require_squarify and not HAS_SQUARIFY:
        print("Aviso: squarify não instalado. Usando visualização alternativa.", file=sys.stderr)
        print("Para treemap completo, instale com: pip install squarify", file=sys.stderr)

    return True


def plot_aa_composition_treemap(aa_composition, sequence_length, output_file='aa_composition_treemap.png'):
    """
    Plota a composição de aminoácidos como um treemap hierárquico usando squarify.

    Args:
        aa_composition: Dicionário {aa: count}
        sequence_length: Comprimento total da sequência
        output_file: Caminho do arquivo de saída
    """
    if not check_dependencies(require_squarify=True):
        return

    # Grupos de aminoácidos por propriedade química
    groups = {
        'Hidrofóbico': ['A', 'V', 'L', 'I', 'M', 'F', 'W', 'P'],
        'Polar': ['S', 'T', 'N', 'Q', 'Y', 'C'],
        'Ácido': ['D', 'E'],
        'Básico': ['K', 'R', 'H'],
        'Glicina': ['G']
    }

    # Cores para cada grupo
    group_colors = {
        'Hidrofóbico': '#E63946',
        'Polar': '#2A9D8F',
        'Ácido': '#F77F00',
        'Básico': '#457B9D',
        'Glicina': '#9D4EDD'
    }

    def get_group(aa):
        """Retorna o grupo do aminoácido."""
        for group, aas in groups.items():
            if aa in aas:
                return group
        return 'Outros'

    # Filtra apenas aminoácidos presentes
    present_aa = {aa: count for aa, count in aa_composition.items() if count > 0}

    if not present_aa:
        print("Erro: Nenhum aminoácido encontrado na sequência.", file=sys.stderr)
        return

    # Ordena por contagem decrescente
    sorted_aa = sorted(present_aa.items(), key=lambda x: -x[1])

    # Prepara dados para o treemap
    sizes = [count for aa, count in sorted_aa]
    labels = [f"{aa}\n{count}\n({count/sequence_length*100:.1f}%)"
              for aa, count in sorted_aa]
    colors = [group_colors.get(get_group(aa), '#888888') for aa, count in sorted_aa]

    # Cria o gráfico
    fig, ax = plt.subplots(figsize=(14, 10))

    if HAS_SQUARIFY:
        # Usa squarify para treemap com quadrados
        squarify.plot(sizes=sizes, label=labels, color=colors,
                     alpha=0.8, ax=ax, text_kwargs={'fontsize': 9, 'weight': 'bold'},
                     edgecolor='white', linewidth=2)
    else:
        # Fallback: gráfico de barras hierárquico
        y_pos = range(len(sorted_aa))
        ax.barh(y_pos, sizes, color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=10)
        ax.set_xlabel('Contagem', fontsize=12, fontweight='bold')

    ax.axis('off')

    # Título
    plt.suptitle('Composição de Aminoácidos - Treemap', fontsize=16, fontweight='bold', y=0.98)

    # Legenda de grupos - PADRONIZADA no lado direito
    legend_elements = [mpatches.Patch(facecolor=color, edgecolor='black', label=group)
                      for group, color in group_colors.items()]
    plt.legend(handles=legend_elements, loc='upper left',
              bbox_to_anchor=(LEGEND_X_POSITION, LEGEND_TOP_Y),
              title='Grupos', **LEGEND_STYLE)

    # Informações estatísticas - PADRONIZADA abaixo da legenda
    total_hydrophobic = sum(aa_composition.get(aa, 0) for aa in groups['Hidrofóbico'])
    total_polar = sum(aa_composition.get(aa, 0) for aa in groups['Polar'])
    total_charged = sum(aa_composition.get(aa, 0) for aa in groups['Ácido'] + groups['Básico'])

    textstr = f'Total: {sequence_length} resíduos\n'
    textstr += f'Hidrofóbicos: {total_hydrophobic/sequence_length*100:.1f}%\n'
    textstr += f'Polares: {total_polar/sequence_length*100:.1f}%\n'
    textstr += f'Carregados: {total_charged/sequence_length*100:.1f}%'

    # Caixa de stats alinhada no lado direito
    ax.text(LEGEND_X_POSITION, INFO_BOX_1_Y, textstr,
           transform=ax.transAxes, fontsize=INFO_BOX_FONT_SIZE,
           verticalalignment='top', bbox=INFO_BOX_STYLE,
           family=INFO_BOX_FAMILY)

    plt.tight_layout(rect=LAYOUT_MARGINS_WITH_LEGEND)
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Treemap de composição salvo em '{output_file}'", file=sys.stderr)
    plt.close()


def plot_aa_composition_bar(aa_composition, sequence_length, output_file='aa_composition_bar.png'):
    """
    Plota a composição de aminoácidos como gráfico de barras.

    Args:
        aa_composition: Dicionário {aa: count}
        sequence_length: Comprimento total da sequência
        output_file: Caminho do arquivo de saída
    """
    if not check_dependencies():
        return

    # Cores por propriedade
    aa_colors = {
        'A': '#E63946', 'V': '#E63946', 'I': '#E63946', 'L': '#E63946',
        'M': '#E63946', 'F': '#E63946', 'W': '#E63946', 'P': '#E63946',
        'S': '#2A9D8F', 'T': '#2A9D8F', 'N': '#2A9D8F', 'Q': '#2A9D8F',
        'Y': '#2A9D8F', 'C': '#2A9D8F',
        'K': '#457B9D', 'R': '#457B9D', 'H': '#457B9D',
        'D': '#F77F00', 'E': '#F77F00',
        'G': '#9D4EDD'
    }

    # Filtra e ordena
    present_aa = {aa: count for aa, count in aa_composition.items() if count > 0}
    sorted_aa = sorted(present_aa.items(), key=lambda x: x[1], reverse=True)

    aa_labels = [aa for aa, _ in sorted_aa]
    aa_percentages = [(count/sequence_length)*100 for _, count in sorted_aa]
    colors = [aa_colors.get(aa, '#888888') for aa in aa_labels]

    # Cria o gráfico
    fig, ax = plt.subplots(figsize=(14, 6))

    bars = ax.bar(range(len(aa_labels)), aa_percentages, color=colors,
                   edgecolor='black', linewidth=0.5, alpha=0.8)

    # Adiciona valores nas barras
    for i, (bar, pct) in enumerate(zip(bars, aa_percentages)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{pct:.1f}%',
                ha='center', va='bottom', fontsize=9, fontweight='bold')

    # Configurações
    ax.set_xlabel('Aminoácido', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequência (%)', fontsize=12, fontweight='bold')
    ax.set_title('Composição de Aminoácidos', fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(aa_labels)))
    ax.set_xticklabels(aa_labels, fontsize=11, fontweight='bold')
    ax.grid(axis='y', alpha=0.3, linestyle=':')
    ax.set_ylim(0, max(aa_percentages) * 1.15)

    # Legenda - PADRONIZADA no lado direito
    legend_elements = [
        mpatches.Patch(facecolor='#E63946', label='Hidrofóbico', edgecolor='black'),
        mpatches.Patch(facecolor='#2A9D8F', label='Polar', edgecolor='black'),
        mpatches.Patch(facecolor='#457B9D', label='Básico', edgecolor='black'),
        mpatches.Patch(facecolor='#F77F00', label='Ácido', edgecolor='black'),
        mpatches.Patch(facecolor='#9D4EDD', label='Glicina', edgecolor='black')
    ]
    ax.legend(handles=legend_elements, loc='upper left',
             bbox_to_anchor=(LEGEND_X_POSITION, LEGEND_TOP_Y),
             title='Grupos', **LEGEND_STYLE)

    plt.tight_layout(rect=LAYOUT_MARGINS_WITH_LEGEND)
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Gráfico de composição salvo em '{output_file}'", file=sys.stderr)
    plt.close()


def plot_hydropathy_profile(sequence, kyte_doolittle_scale, output_file='hydropathy_profile.png', window=9):
    """
    Plota o perfil de hidrofobicidade usando janela deslizante (Kyte-Doolittle).

    Args:
        sequence: Sequência de aminoácidos
        kyte_doolittle_scale: Dicionário com escala de Kyte-Doolittle
        output_file: Caminho do arquivo de saída
        window: Tamanho da janela deslizante (padrão: 9)
    """
    if not check_dependencies():
        return

    if len(sequence) < window:
        print(f"Aviso: Sequência muito curta ({len(sequence)} aa) para janela de {window}.", file=sys.stderr)
        window = max(1, len(sequence) // 2)

    # Calcula hidrofobicidade em janela deslizante
    scores = []
    positions = []

    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i+window]
        score = sum(kyte_doolittle_scale.get(aa, 0) for aa in window_seq) / window
        scores.append(score)
        positions.append(i + window//2 + 1)

    # Cria o gráfico
    fig, ax = plt.subplots(figsize=(14, 6))

    # Plot da linha
    ax.plot(positions, scores, linewidth=2.5, color='#264653', label='Hidrofobicidade', zorder=3)
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5, zorder=1)

    # Preenche áreas hidrofóbicas e hidrofílicas
    ax.fill_between(positions, scores, 0, where=[s > 0 for s in scores],
                     alpha=0.4, color='#E63946', label='Hidrofóbico', interpolate=True, zorder=2)
    ax.fill_between(positions, scores, 0, where=[s < 0 for s in scores],
                     alpha=0.4, color='#2A9D8F', label='Hidrofílico', interpolate=True, zorder=2)

    # Configurações
    ax.set_xlabel('Posição na Sequência', fontsize=12, fontweight='bold')
    ax.set_ylabel('Hidrofobicidade (Kyte-Doolittle)', fontsize=12, fontweight='bold')
    ax.set_title(f'Perfil de Hidrofobicidade (janela = {window} resíduos)', fontsize=14, fontweight='bold')

    # Estatísticas - calculadas primeiro
    avg_score = sum(scores) / len(scores)
    max_score = max(scores)
    min_score = min(scores)
    max_pos = positions[scores.index(max_score)]
    min_pos = positions[scores.index(min_score)]

    # Legenda - PADRONIZADA no lado direito
    legend_items = ax.legend(loc='upper left',
                            bbox_to_anchor=(LEGEND_X_POSITION, LEGEND_TOP_Y),
                            **LEGEND_STYLE)

    # Nota sobre regiões transmembrana - PADRONIZADA abaixo da legenda
    if max_score > 1.6:
        note_text = "Nota:\nPicos > 1.6 sugerem\npossíveis domínios\ntransmembrana"
        ax.text(LEGEND_X_POSITION, INFO_BOX_1_Y, note_text,
               transform=ax.transAxes, fontsize=INFO_BOX_FONT_SIZE,
               verticalalignment='top', bbox=INFO_BOX_STYLE,
               family=INFO_BOX_FAMILY, linespacing=1.5)

    # Estatísticas - PADRONIZADA mais abaixo no lado direito
    textstr = f'Média: {avg_score:.3f}\n'
    textstr += f'Máx: {max_score:.3f} (pos {max_pos})\n'
    textstr += f'Mín: {min_score:.3f} (pos {min_pos})'

    ax.text(LEGEND_X_POSITION, INFO_BOX_2_Y, textstr,
           transform=ax.transAxes, fontsize=INFO_BOX_FONT_SIZE,
           verticalalignment='top', bbox=INFO_BOX_STYLE,
           family=INFO_BOX_FAMILY)

    ax.grid(True, alpha=0.3, linestyle=':', zorder=0)
    ax.set_xlim(positions[0], positions[-1])

    plt.tight_layout(rect=LAYOUT_MARGINS_WITH_LEGEND)
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Perfil de hidrofobicidade salvo em '{output_file}'", file=sys.stderr)
    plt.close()


def plot_sasa_profile(sasa_per_atom, output_file='sasa_profile.png'):
    """
    Plota o perfil de SASA médio por resíduo ao longo da sequência.

    Args:
        sasa_per_atom: Lista de dicionários com chain_id, res_num, sasa
        output_file: Caminho do arquivo de saída
    """
    if not check_dependencies(require_numpy=True):
        return

    from collections import defaultdict

    # Calcula SASA médio por resíduo
    residue_sasa = defaultdict(list)
    for atom in sasa_per_atom:
        residue_key = (atom['chain_id'], atom['res_num'])
        residue_sasa[residue_key].append(atom['sasa'])

    # Ordena por número de resíduo
    sorted_residues = sorted(residue_sasa.keys(), key=lambda x: x[1])
    residue_numbers = [res_num for _, res_num in sorted_residues]
    avg_sasa = [sum(residue_sasa[key])/len(residue_sasa[key]) for key in sorted_residues]

    # Cria o gráfico
    fig, ax = plt.subplots(figsize=(14, 6))

    # Plot com cores baseadas em SASA
    scatter = ax.scatter(residue_numbers, avg_sasa, c=avg_sasa, cmap='RdYlBu_r',
                         s=40, alpha=0.7, edgecolors='black', linewidths=0.5, zorder=3)
    ax.plot(residue_numbers, avg_sasa, linewidth=1.5, color='#264653', alpha=0.6, zorder=2)

    # Threshold para exposto vs enterrado (20 Ų é típico)
    threshold = 20.0
    ax.axhline(y=threshold, color='red', linestyle='--', linewidth=1.5,
               label=f'Threshold ({threshold:.0f} Ų)', alpha=0.7, zorder=1)

    # Preenche áreas
    ax.fill_between(residue_numbers, avg_sasa, threshold,
                     where=[s > threshold for s in avg_sasa],
                     alpha=0.2, color='#2A9D8F', label='Exposto', zorder=0)
    ax.fill_between(residue_numbers, avg_sasa, threshold,
                     where=[s <= threshold for s in avg_sasa],
                     alpha=0.2, color='#E76F51', label='Enterrado', zorder=0)

    # Configurações
    ax.set_xlabel('Número do Resíduo', fontsize=12, fontweight='bold')
    ax.set_ylabel('SASA Médio (Ų)', fontsize=12, fontweight='bold')
    ax.set_title('Perfil de Acessibilidade ao Solvente (SASA)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle=':')

    # Colorbar - ANTES da legenda para posicionar corretamente
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label('SASA (Ų)', fontsize=10, fontweight='bold')

    # Legenda - Deslocada para direita para não sobrepor o colorbar
    ax.legend(loc='upper left',
             bbox_to_anchor=(LEGEND_X_WITH_COLORBAR, LEGEND_TOP_Y),
             **LEGEND_STYLE)

    # Estatísticas - Deslocadas para direita, abaixo da legenda
    num_exposed = sum(1 for s in avg_sasa if s > threshold)
    num_buried = len(avg_sasa) - num_exposed
    mean_sasa = sum(avg_sasa) / len(avg_sasa)
    max_sasa = max(avg_sasa)
    max_pos = residue_numbers[avg_sasa.index(max_sasa)]

    textstr = f'Total: {len(avg_sasa)} resíduos\n'
    textstr += f'Expostos (>{threshold:.0f} Ų): {num_exposed}\n'
    textstr += f'Enterrados (≤{threshold:.0f} Ų): {num_buried}\n'
    textstr += f'SASA médio: {mean_sasa:.2f} Ų\n'
    textstr += f'SASA máx: {max_sasa:.2f} Ų (res {max_pos})'

    ax.text(LEGEND_X_WITH_COLORBAR, INFO_BOX_1_Y, textstr,
           transform=ax.transAxes, fontsize=INFO_BOX_FONT_SIZE,
           verticalalignment='top', bbox=INFO_BOX_STYLE,
           family=INFO_BOX_FAMILY)

    plt.tight_layout(rect=LAYOUT_MARGINS_WITH_COLORBAR)
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Perfil de SASA salvo em '{output_file}'", file=sys.stderr)
    plt.close()


def plot_contact_map(contacts, max_residue, output_file='contact_map.png', threshold=8.0):
    """
    Plota um mapa de contatos intramoleculares.

    Args:
        contacts: Lista de tuplas (res1, res2, distance)
        max_residue: Número máximo de resíduo
        output_file: Caminho do arquivo de saída
        threshold: Distância threshold usada
    """
    if not check_dependencies(require_numpy=True):
        return

    # Cria matriz de distâncias (não binária)
    contact_matrix = np.zeros((max_residue, max_residue))

    # Lista para rastrear min/max das distâncias
    distances = []

    for res1, res2, dist in contacts:
        i, j = res1 - 1, res2 - 1
        if 0 <= i < max_residue and 0 <= j < max_residue:
            contact_matrix[i, j] = dist
            contact_matrix[j, i] = dist
            distances.append(dist)

    # Calcula min/max das distâncias para a escala
    if distances:
        min_dist = min(distances)
        max_dist = max(distances)
    else:
        min_dist, max_dist = 0, threshold

    # Cria o gráfico
    fig, ax = plt.subplots(figsize=(10, 10))

    # Heatmap com escala de distâncias
    # Usar máscara para não mostrar zeros (sem contato) como azul escuro
    masked_matrix = np.ma.masked_where(contact_matrix == 0, contact_matrix)

    im = ax.imshow(masked_matrix, cmap='viridis_r', origin='lower',
                   interpolation='nearest', aspect='auto',
                   vmin=min_dist, vmax=max_dist)

    # Linha diagonal
    ax.plot([0, max_residue-1], [0, max_residue-1], 'r--', linewidth=1, alpha=0.5, label='Diagonal')

    # Configurações
    ax.set_xlabel('Resíduo i', fontsize=12, fontweight='bold')
    ax.set_ylabel('Resíduo j', fontsize=12, fontweight='bold')
    ax.set_title(f'Mapa de Contatos Intramoleculares\n({len(contacts)} contatos)',
                 fontsize=14, fontweight='bold')

    # Colorbar com escala de distâncias
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Distância (Å)', fontsize=10, fontweight='bold')

    # Define ticks mostrando min, max e valores intermediários
    if distances:
        tick_values = np.linspace(min_dist, max_dist, 5)
        cbar.set_ticks(tick_values)
        cbar.set_ticklabels([f'{val:.1f}' for val in tick_values])

    # Grade
    ax.grid(True, which='both', alpha=0.2, linestyle=':', color='gray')

    # Estatísticas - PADRONIZADA no lado direito
    total_possible = (max_residue * (max_residue - 1)) // 2
    contact_density = (len(contacts) / total_possible) * 100 if total_possible > 0 else 0

    textstr = f'Resíduos totais: {max_residue}\n'
    textstr += f'Contatos: {len(contacts)}\n'
    textstr += f'Densidade: {contact_density:.2f}%\n'
    textstr += f'Threshold: {threshold:.1f} Å'

    ax.text(LEGEND_X_POSITION, LEGEND_TOP_Y, textstr,
           transform=ax.transAxes, fontsize=INFO_BOX_FONT_SIZE,
           verticalalignment='top', bbox=INFO_BOX_STYLE,
           family=INFO_BOX_FAMILY)

    plt.tight_layout(rect=LAYOUT_MARGINS_WITH_LEGEND)
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Mapa de contatos salvo em '{output_file}'", file=sys.stderr)
    plt.close()


if __name__ == "__main__":
    print("Este módulo contém funções de visualização para o BioHub.")
    print("Não deve ser executado diretamente.")
    print("\nDependências necessárias:")
    print("  - matplotlib: ", "OK" if HAS_MATPLOTLIB else "AUSENTE")
    print("  - numpy: ", "OK" if HAS_NUMPY else "AUSENTE")
    print("  - squarify: ", "OK" if HAS_SQUARIFY else "AUSENTE (opcional)")
