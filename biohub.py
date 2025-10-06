#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Uma ferramenta de linha de comando para análises bioinformáticas básicas,
com o mínimo de dependências externas.
"""

import sys
import math
import argparse

# --- Constantes e Dicionários de Dados ---
# Mapeamento de código de 3 letras para 1 letra de aminoácidos
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# Pesos moleculares (média, em Daltons) para cada resíduo
MOLECULAR_WEIGHT = {
    'A': 89.09,  'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16,
    'E': 147.13, 'Q': 146.15, 'G': 75.07,  'H': 155.16, 'I': 131.17,
    'L': 131.17, 'K': 146.19, 'M': 149.21, 'F': 165.19, 'P': 115.13,
    'S': 105.09, 'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15
}

# Escala de hidropatia de Kyte & Doolittle
KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Valores de pKa para grupos ionizáveis
PKA_VALUES = {
    'C-term': 3.65,
    'N-term': 8.0,
    'D': 3.9, 'E': 4.07, 'H': 6.5, 'C': 8.5,
    'Y': 10.0, 'K': 10.0, 'R': 12.0
}

# --- Funções Principais ---

def get_sequence_from_pdb(pdb_filepath: str) -> str:
    """
    Extrai a sequência de aminoácidos de um arquivo PDB.
    Considera apenas os registros ATOM da primeira cadeia (chain) encontrada.
    
    Args:
        pdb_filepath: O caminho para o arquivo PDB.

    Returns:
        A sequência de aminoácidos em formato de uma letra.
    """
    sequence = ""
    processed_residues = set()
    first_chain_id = None

    try:
        with open(pdb_filepath, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    current_chain_id = line[21]
                    # Define a cadeia a ser processada com base na primeira encontrada
                    if first_chain_id is None:
                        first_chain_id = current_chain_id
                    
                    # Processa apenas a primeira cadeia
                    if current_chain_id == first_chain_id:
                        residue_num = int(line[22:26])
                        residue_name = line[17:20]

                        # Evita duplicatas para o mesmo resíduo
                        if residue_num not in processed_residues:
                            if residue_name in THREE_TO_ONE:
                                sequence += THREE_TO_ONE[residue_name]
                                processed_residues.add(residue_num)
    except FileNotFoundError:
        print(f"Erro: Arquivo não encontrado em '{pdb_filepath}'", file=sys.stderr)
        return ""
    except Exception as e:
        print(f"Erro ao processar o arquivo PDB: {e}", file=sys.stderr)
        return ""
        
    return sequence

def calculate_physicochemical_properties(sequence: str):
    """
    Calcula e exibe as propriedades físico-químicas de uma sequência de aminoácidos.
    
    Args:
        sequence: A sequência de aminoácidos.
    """
    if not sequence:
        print("Erro: A sequência está vazia.", file=sys.stderr)
        return

    # 1. Composição de aminoácidos
    aa_composition = {aa: sequence.count(aa) for aa in MOLECULAR_WEIGHT.keys()}
    
    # 2. Peso Molecular
    # Soma dos pesos dos resíduos - (N-1) * peso da água
    mw = sum(MOLECULAR_WEIGHT[aa] * aa_composition[aa] for aa in aa_composition)
    water_weight = 18.015
    mw -= (len(sequence) - 1) * water_weight

    # 3. Ponto Isoelétrico (pI)
    pi = 0.0
    min_charge_diff = float('inf')
    
    # Itera sobre uma faixa de pH para encontrar onde a carga líquida é próxima de zero
    for ph in [i * 0.01 for i in range(1401)]:
        # Carga do N-terminal (positivo)
        net_charge = (10**PKA_VALUES['N-term']) / (10**PKA_VALUES['N-term'] + 10**ph)
        # Carga do C-terminal (negativo)
        net_charge -= (10**ph) / (10**PKA_VALUES['C-term'] + 10**ph)
        
        # Cargas dos resíduos
        for aa, count in aa_composition.items():
            if aa in ['R', 'H', 'K']: # Resíduos básicos (carga positiva)
                net_charge += count * (10**PKA_VALUES[aa]) / (10**PKA_VALUES[aa] + 10**ph)
            elif aa in ['D', 'E', 'C', 'Y']: # Resíduos ácidos (carga negativa)
                net_charge -= count * (10**ph) / (10**PKA_VALUES[aa] + 10**ph)

        if abs(net_charge) < min_charge_diff:
            min_charge_diff = abs(net_charge)
            pi = ph

    # 4. Coeficiente de extinção
    n_Y = aa_composition['Y']
    n_W = aa_composition['W']
    n_C = aa_composition['C']
    ext_coeff_reduced = (n_W * 5500) + (n_Y * 1490)
    ext_coeff_oxidized = ext_coeff_reduced + (n_C // 2) * 125

    # 5. GRAVY (Grand Average of Hydropathicity)
    gravy = sum(KYTE_DOOLITTLE.get(aa, 0) for aa in sequence) / len(sequence)

    print("--- Propriedades Físico-Químicas ---")
    print(f"Comprimento da Sequência: {len(sequence)}")
    print(f"Peso Molecular (MW): {mw:.2f} Da")
    print(f"Ponto Isoelétrico (pI): {pi:.2f}")
    print(f"Coeficiente de Extinção (reduzido): {ext_coeff_reduced} M⁻¹cm⁻¹")
    print(f"Coeficiente de Extinção (oxidado): {ext_coeff_oxidized} M⁻¹cm⁻¹")
    print(f"GRAVY: {gravy:.3f}")
    print("\nComposição de Aminoácidos:")
    for aa, count in sorted(aa_composition.items()):
        if count > 0:
            print(f"  {aa}: {count} ({(count / len(sequence)) * 100:.2f}%)")


def calculate_intramolecular_contacts(pdb_filepath: str, threshold: float = 8.0):
    """
    Calcula e lista os contatos intramoleculares com base na distância
    entre os átomos de Carbono Alfa (CA).
    
    Args:
        pdb_filepath: O caminho para o arquivo PDB.
        threshold: A distância máxima (em Angstroms) para considerar um contato.
    """
    ca_atoms = {}
    first_chain_id = None

    try:
        with open(pdb_filepath, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM") and line[13:16].strip() == "CA":
                    current_chain_id = line[21]
                    if first_chain_id is None:
                        first_chain_id = current_chain_id
                    
                    if current_chain_id == first_chain_id:
                        res_num = int(line[22:26])
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ca_atoms[res_num] = (x, y, z)
    except FileNotFoundError:
        print(f"Erro: Arquivo não encontrado em '{pdb_filepath}'", file=sys.stderr)
        return

    if not ca_atoms:
        print("Nenhum átomo de Carbono Alfa encontrado na primeira cadeia do arquivo.", file=sys.stderr)
        return

    residues = sorted(ca_atoms.keys())
    contacts = []

    # Itera sobre todos os pares de resíduos sem repetição
    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            res1_num = residues[i]
            res2_num = residues[j]
            
            # Ignora resíduos adjacentes
            if abs(res1_num - res2_num) <= 1:
                continue

            coord1 = ca_atoms[res1_num]
            coord2 = ca_atoms[res2_num]
            
            # Calcula a distância euclidiana
            distance = math.sqrt(
                (coord1[0] - coord2[0])**2 +
                (coord1[1] - coord2[1])**2 +
                (coord1[2] - coord2[2])**2
            )

            if distance <= threshold:
                contacts.append((res1_num, res2_num, distance))

    print(f"--- Contatos Intramoleculares (Limiar = {threshold:.1f} Å) ---")
    if not contacts:
        print("Nenhum contato encontrado com os critérios fornecidos.")
    else:
        print(f"Total de contatos encontrados: {len(contacts)}")
        print("Resíduo 1 | Resíduo 2 | Distância (Å)")
        print("----------|-----------|---------------")
        for c in contacts:
            print(f"{c[0]:<9} | {c[1]:<9} | {c[2]:.3f}")

def predict_solvent_exposure(pdb_filepath: str, window_size: int = 9):
    """
    Prevê regiões expostas ao solvente usando o método de Kyte-Doolittle
    com uma janela deslizante.
    
    Args:
        pdb_filepath: O caminho para o arquivo PDB.
        window_size: O tamanho da janela deslizante para calcular a média de hidropatia.
    """
    if window_size % 2 == 0:
        print("Erro: O tamanho da janela deve ser um número ímpar.", file=sys.stderr)
        return
        
    sequence = get_sequence_from_pdb(pdb_filepath)
    if not sequence:
        return
        
    scores = []
    half_window = window_size // 2

    # Itera sobre a sequência para calcular a pontuação de cada resíduo
    for i in range(len(sequence)):
        start = max(0, i - half_window)
        end = min(len(sequence), i + half_window + 1)
        window_seq = sequence[start:end]
        
        # Calcula a média da hidropatia na janela
        window_score = sum(KYTE_DOOLITTLE.get(aa, 0) for aa in window_seq) / len(window_seq)
        scores.append(window_score)

    print(f"--- Predição de Exposição ao Solvente (Kyte-Doolittle, Janela={window_size}) ---")
    print("Scores > 0 indicam regiões hidrofóbicas (provavelmente internas).")
    print("Scores < 0 indicam regiões hidrofílicas (provavelmente expostas).")
    print("\nPos. | Res. | Score de Hidropatia")
    print("-----|------|---------------------")
    for i, score in enumerate(scores):
        print(f"{i+1:<4} | {sequence[i]:<4} | {score:.3f}")

# --- Configuração da Interface de Linha de Comando ---

def main():
    parser = argparse.ArgumentParser(
        description="Ferramenta de bioinformática para análise de PDB e sequências."
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Função a ser executada")

    # Subcomando para converter PDB para FASTA
    parser_fasta = subparsers.add_parser("fasta", help="Converte um arquivo PDB em uma sequência FASTA.")
    parser_fasta.add_argument("pdb_file", help="Caminho para o arquivo PDB de entrada.")
    parser_fasta.add_argument("-o", "--output", help="Caminho para o arquivo FASTA de saída (opcional, padrão: stdout).")

    # Subcomando para propriedades físico-químicas
    parser_physchem = subparsers.add_parser("physchem", help="Calcula propriedades físico-químicas de uma sequência.")
    parser_physchem.add_argument("sequence", help="A sequência de aminoácidos a ser analisada.")
    
    # Subcomando para contatos intramoleculares
    parser_contacts = subparsers.add_parser("contacts", help="Calcula contatos intramoleculares a partir de um arquivo PDB.")
    parser_contacts.add_argument("pdb_file", help="Caminho para o arquivo PDB de entrada.")
    parser_contacts.add_argument("-t", "--threshold", type=float, default=8.0, help="Distância máxima em Angstroms para um contato (padrão: 8.0).")

    # Subcomando para predição de exposição ao solvente
    parser_exposure = subparsers.add_parser("exposure", help="Prevê a exposição ao solvente com a escala Kyte-Doolittle.")
    parser_exposure.add_argument("pdb_file", help="Caminho para o arquivo PDB de entrada.")
    parser_exposure.add_argument("-w", "--window", type=int, default=9, help="Tamanho da janela deslizante (padrão: 9).")
    
    args = parser.parse_args()

    # Executa a função correspondente ao comando
    if args.command == "fasta":
        sequence = get_sequence_from_pdb(args.pdb_file)
        if sequence:
            header = f">sequence_from_{args.pdb_file.split('/')[-1]}"
            fasta_output = f"{header}\n{sequence}\n"
            if args.output:
                with open(args.output, 'w') as f:
                    f.write(fasta_output)
                print(f"Sequência FASTA salva em '{args.output}'")
            else:
                print(fasta_output)
    elif args.command == "physchem":
        calculate_physicochemical_properties(args.sequence.upper())
    elif args.command == "contacts":
        calculate_intramolecular_contacts(args.pdb_file, args.threshold)
    elif args.command == "exposure":
        predict_solvent_exposure(args.pdb_file, args.window)

if __name__ == "__main__":
    main()