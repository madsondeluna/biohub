#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioHub: Uma ferramenta de linha de comando para análises bioinformáticas,
com o mínimo de dependências externas.
"""

import sys
import math
import argparse
import subprocess
import os
import tempfile

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
    'C-term': 3.65, 'N-term': 8.0, 'D': 3.9, 'E': 4.07, 'H': 6.5, 
    'C': 8.5, 'Y': 10.0, 'K': 10.0, 'R': 12.0
}

# Raios de Van der Waals (em Angstroms) para tipos de átomos comuns
# Extraído da força-campo AMBER
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'S': 1.80, 'P': 1.80, 'F': 1.47, 'CL': 1.75,
    'BR': 1.85, 'I': 1.98, 'DEFAULT': 1.70
}

# --- Funções Auxiliares ---

def parse_pdb_atoms(pdb_filepath: str):
    """Lê um arquivo PDB e retorna uma lista de átomos com suas propriedades."""
    atoms = []
    try:
        with open(pdb_filepath, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom = {
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                        "res_name": line[17:20].strip(),
                        "res_num": int(line[22:26]),
                        "chain_id": line[21],
                        "element": line[76:78].strip().upper() or line[12:14].strip().upper()
                    }
                    atoms.append(atom)
    except FileNotFoundError:
        print(f"Erro: Arquivo não encontrado em '{pdb_filepath}'", file=sys.stderr)
        return []
    return atoms

def generate_sphere_points(n_points: int):
    """Gera pontos uniformemente distribuídos na superfície de uma esfera unitária."""
    points = []
    phi = (1 + math.sqrt(5)) / 2  # Golden ratio
    for i in range(n_points):
        y = 1 - (2 * i / (n_points - 1))
        radius = math.sqrt(1 - y * y)
        theta = 2 * math.pi * i / phi
        x = math.cos(theta) * radius
        z = math.sin(theta) * radius
        points.append((x, y, z))
    return points


# --- Funções Principais ---

def get_sequence_from_pdb(pdb_filepath: str) -> str:
    """Extrai a sequência de aminoácidos de um arquivo PDB."""
    sequence = ""
    processed_residues = set()
    first_chain_id = None
    
    try:
        with open(pdb_filepath, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    current_chain_id = line[21]
                    if first_chain_id is None:
                        first_chain_id = current_chain_id
                    
                    if current_chain_id == first_chain_id:
                        res_id = (int(line[22:26]), line[21])
                        if res_id not in processed_residues:
                            res_name = line[17:20]
                            if res_name in THREE_to_ONE:
                                sequence += THREE_TO_ONE[res_name]
                                processed_residues.add(res_id)
    except FileNotFoundError:
        print(f"Erro: Arquivo não encontrado em '{pdb_filepath}'", file=sys.stderr)
        return ""
    return sequence

def calculate_physicochemical_properties(sequence: str):
    """Calcula e exibe as propriedades físico-químicas de uma sequência."""
    # (Implementação mantida da versão anterior)
    if not sequence:
        print("Erro: A sequência está vazia.", file=sys.stderr)
        return
    aa_composition = {aa: sequence.count(aa) for aa in MOLECULAR_WEIGHT.keys()}
    mw = sum(MOLECULAR_WEIGHT[aa] * aa_composition[aa] for aa in aa_composition)
    mw -= (len(sequence) - 1) * 18.015
    pi = 0.0
    min_charge_diff = float('inf')
    for ph_int in range(1401):
        ph = ph_int * 0.01
        net_charge = (10**PKA_VALUES['N-term']) / (10**PKA_VALUES['N-term'] + 10**ph)
        net_charge -= (10**ph) / (10**PKA_VALUES['C-term'] + 10**ph)
        for aa, count in aa_composition.items():
            if aa in ['R', 'H', 'K']:
                net_charge += count * (10**PKA_VALUES[aa]) / (10**PKA_VALUES[aa] + 10**ph)
            elif aa in ['D', 'E', 'C', 'Y']:
                net_charge -= count * (10**ph) / (10**PKA_VALUES[aa] + 10**ph)
        if abs(net_charge) < min_charge_diff:
            min_charge_diff = abs(net_charge)
            pi = ph
    n_Y, n_W, n_C = aa_composition['Y'], aa_composition['W'], aa_composition['C']
    ext_coeff_reduced = (n_W * 5500) + (n_Y * 1490)
    ext_coeff_oxidized = ext_coeff_reduced + (n_C // 2) * 125
    gravy = sum(KYTE_DOOLITTLE.get(aa, 0) for aa in sequence) / len(sequence)

    print("--- Propriedades Físico-Químicas ---")
    print(f"Comprimento da Sequência: {len(sequence)}")
    print(f"Peso Molecular (MW): {mw:.2f} Da")
    print(f"Ponto Isoelétrico (pI): {pi:.2f}")
    print(f"Coeficiente de Extinção (reduzido): {ext_coeff_reduced} M⁻¹cm⁻¹")
    print(f"Coeficiente de Extinção (oxidado): {ext_coeff_oxidized} M⁻¹cm⁻¹")
    print(f"GRAVY: {gravy:.3f}")

def calculate_intramolecular_contacts(pdb_filepath: str, threshold: float = 8.0):
    """Calcula contatos intramoleculares com base na distância entre Carbonos Alfa."""
    # (Implementação mantida da versão anterior)
    ca_atoms = {}
    atoms = parse_pdb_atoms(pdb_filepath)
    if not atoms: return
    
    first_chain_id = atoms[0]["chain_id"]
    with open(pdb_filepath, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[13:16].strip() == "CA" and line[21] == first_chain_id:
                 ca_atoms[int(line[22:26])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

    residues = sorted(ca_atoms.keys())
    contacts = []
    for i in range(len(residues)):
        for j in range(i + 1, len(residues)):
            res1_num, res2_num = residues[i], residues[j]
            if abs(res1_num - res2_num) <= 1: continue
            coord1, coord2 = ca_atoms[res1_num], ca_atoms[res2_num]
            distance = math.sqrt(sum((c1 - c2)**2 for c1, c2 in zip(coord1, coord2)))
            if distance <= threshold:
                contacts.append((res1_num, res2_num, distance))

    print(f"--- Contatos Intramoleculares (Limiar = {threshold:.1f} Å) ---")
    if contacts:
        print(f"Total de contatos: {len(contacts)}")
        for c in contacts: print(f"Res {c[0]} - Res {c[1]}: {c[2]:.3f} Å")
    else:
        print("Nenhum contato encontrado.")

def predict_solvent_exposure(pdb_filepath: str, window_size: int = 9):
    """Prevê a exposição ao solvente usando a escala Kyte-Doolittle."""
    # (Implementação mantida da versão anterior)
    if window_size % 2 == 0:
        print("Erro: O tamanho da janela deve ser ímpar.", file=sys.stderr)
        return
    sequence = get_sequence_from_pdb(pdb_filepath)
    if not sequence: return
    
    scores = []
    half_window = window_size // 2
    for i in range(len(sequence)):
        window_seq = sequence[max(0, i - half_window) : min(len(sequence), i + half_window + 1)]
        scores.append(sum(KYTE_DOOLITTLE.get(aa, 0) for aa in window_seq) / len(window_seq))
    
    print(f"--- Predição de Exposição (Kyte-Doolittle, Janela={window_size}) ---")
    print("Pos. | Res. | Score de Hidropatia")
    for i, score in enumerate(scores):
        print(f"{i+1:<4} | {sequence[i]:<4} | {score:.3f}")

def calculate_sasa(pdb_filepath: str, probe_radius: float, n_points: int):
    """Calcula a Área de Superfície Acessível ao Solvente (SASA) via algoritmo de Shrake-Rupley."""
    atoms = parse_pdb_atoms(pdb_filepath)
    if not atoms: return
    
    sphere_points = generate_sphere_points(n_points)
    
    sasa_per_residue = {}
    total_sasa = 0.0

    for i, atom_i in enumerate(atoms):
        radius_i = VDW_RADII.get(atom_i["element"], VDW_RADII['DEFAULT'])
        extended_radius = radius_i + probe_radius
        
        accessible_points = 0
        for sp in sphere_points:
            point_is_accessible = True
            point = (atom_i["x"] + extended_radius * sp[0],
                     atom_i["y"] + extended_radius * sp[1],
                     atom_i["z"] + extended_radius * sp[2])
            
            for j, atom_j in enumerate(atoms):
                if i == j: continue
                
                radius_j = VDW_RADII.get(atom_j["element"], VDW_RADII['DEFAULT'])
                dist_sq = sum((p - c)**2 for p, c in zip(point, (atom_j["x"], atom_j["y"], atom_j["z"])))
                
                if dist_sq < (radius_j + probe_radius)**2:
                    point_is_accessible = False
                    break
            
            if point_is_accessible:
                accessible_points += 1
        
        atom_sasa = (accessible_points / n_points) * 4.0 * math.pi * extended_radius**2
        total_sasa += atom_sasa
        
        res_id = (atom_i["res_num"], atom_i["res_name"])
        sasa_per_residue[res_id] = sasa_per_residue.get(res_id, 0) + atom_sasa

    print(f"--- Análise de SASA (Shrake-Rupley) ---")
    print(f"Raio da sonda: {probe_radius} Å, Pontos por átomo: {n_points}")
    print(f"\nSASA Total da Molécula: {total_sasa:.2f} Å²")
    print("\nSASA por Resíduo:")
    print("Resíduo  | SASA (Å²)")
    print("---------|----------")
    for res_id, sasa in sorted(sasa_per_residue.items()):
        print(f"{res_id[1]:<3} {res_id[0]:<4} | {sasa:.2f}")

def run_apbs_analysis(pdb_filepath: str, no_cleanup: bool):
    """Wrapper para executar PDB2PQR e APBS e calcular a energia de solvatação."""
    print("--- Análise de Energia de Solvatação com APBS ---")
    
    for exe in ["pdb2pqr", "apbs"]:
        if not any(os.access(os.path.join(path, exe), os.X_OK) for path in os.environ.get("PATH", "").split(os.pathsep)):
            print(f"Erro: O executável '{exe}' não foi encontrado no PATH do sistema.", file=sys.stderr)
            print("Por favor, instale PDB2PQR e APBS e garanta que estejam acessíveis.", file=sys.stderr)
            return

    work_dir = tempfile.mkdtemp()
    try:
        base_name = os.path.splitext(os.path.basename(pdb_filepath))[0]
        pqr_path = os.path.join(work_dir, f"{base_name}.pqr")
        
        print(f"1. Executando pdb2pqr em '{pdb_filepath}'...")
        try:
            subprocess.run(
                ["pdb2pqr", "--ff=amber", pdb_filepath, pqr_path],
                check=True, capture_output=True, text=True, timeout=60
            )
        except subprocess.CalledProcessError as e:
            print(f"Erro ao executar pdb2pqr:\n{e.stderr}", file=sys.stderr)
            return
            
        apbs_in_path = os.path.join(work_dir, "apbs.in")
        apbs_input_content = f"""
read
    mol pqr {os.path.basename(pqr_path)}
end
elec
    mg-auto
end
quit
"""
        with open(apbs_in_path, 'w') as f: f.write(apbs_input_content)
        print("2. Arquivo de input do APBS ('apbs.in') criado.")

        print("3. Executando APBS... (Isso pode levar alguns instantes)")
        try:
            result = subprocess.run(
                ["apbs", "apbs.in"],
                check=True, capture_output=True, text=True, cwd=work_dir, timeout=300
            )
        except subprocess.CalledProcessError as e:
            print(f"Erro ao executar APBS:\n{e.stderr}", file=sys.stderr)
            return

        print("4. Analisando os resultados...")
        energy = "Não encontrada"
        for line in result.stdout.splitlines():
            if "Global net ELEC energy" in line:
                try:
                    energy = f"{float(line.split()[-2]):.4f} kJ/mol"
                    break
                except (ValueError, IndexError):
                    energy = "Erro ao analisar valor"
        
        print("\n--- Resultado da Energia de Solvatação Eletrostática ---")
        print(f"Energia: {energy}")

    finally:
        if no_cleanup:
            print(f"\nArquivos intermediários mantidos em: {work_dir}")
        else:
            import shutil
            shutil.rmtree(work_dir)

def main():
    parser = argparse.ArgumentParser(
        description="BioHub: Ferramenta de bioinformática para análise de PDB e sequências.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Função a ser executada")

    parser_fasta = subparsers.add_parser("fasta", help="Converte um arquivo PDB em uma sequência FASTA.")
    parser_fasta.add_argument("pdb_file", help="Caminho para o arquivo PDB de entrada.")
    parser_fasta.add_argument("-o", "--output", help="Arquivo de saída (padrão: stdout).")

    parser_physchem = subparsers.add_parser("physchem", help="Calcula propriedades físico-químicas de uma sequência.")
    parser_physchem.add_argument("sequence", help="A sequência de aminoácidos a ser analisada.")

    parser_contacts = subparsers.add_parser("contacts", help="Calcula contatos intramoleculares a partir de um arquivo PDB.")
    parser_contacts.add_argument("pdb_file", help="Caminho para o arquivo PDB de entrada.")
    parser_contacts.add_argument("-t", "--threshold", type=float, default=8.0, help="Distância máxima em Angstroms para um contato (padrão: 8.0).")

    parser_exposure = subparsers.add_parser("exposure", help="Prevê a exposição ao solvente com a escala Kyte-Doolittle.")
    parser_exposure.add_argument("pdb_file", help="Caminho para o arquivo PDB de entrada.")
    parser_exposure.add_argument("-w", "--window", type=int, default=9, help="Tamanho da janela deslizante (padrão: 9).")

    parser_sasa = subparsers.add_parser("sasa", help="Calcula a Área de Superfície Acessível ao Solvente (SASA).")
    parser_sasa.add_argument("pdb_file", help="Caminho para o arquivo PDB de entrada.")
    parser_sasa.add_argument("--probe-radius", type=float, default=1.4, help="Raio da sonda do solvente em Angstroms (padrão: 1.4 para água).")
    parser_sasa.add_argument("--num-points", type=int, default=960, help="Número de pontos na superfície de cada átomo para o cálculo (padrão: 960).")

    parser_apbs = subparsers.add_parser("apbs", help="Calcula a energia de solvatação eletrostática via APBS (requer PDB2PQR e APBS).")
    parser_apbs.add_argument("pdb_file", help="Caminho para o arquivo PDB de entrada.")
    parser_apbs.add_argument("--no-cleanup", action="store_true", help="Não remove os arquivos intermediários (PQR, apbs.in, etc.) após a execução.")

    args = parser.parse_args()

    command_functions = {
        "fasta": lambda a: handle_fasta(a.pdb_file, a.output),
        "physchem": lambda a: calculate_physicochemical_properties(a.sequence.upper()),
        "contacts": lambda a: calculate_intramolecular_contacts(a.pdb_file, a.threshold),
        "exposure": lambda a: predict_solvent_exposure(a.pdb_file, a.window),
        "sasa": lambda a: calculate_sasa(a.pdb_file, a.probe_radius, a.num_points),
        "apbs": lambda a: run_apbs_analysis(a.pdb_file, a.no_cleanup)
    }
    
    # Executa a função correspondente
    command_functions[args.command](args)

def handle_fasta(pdb_file, output_file):
    sequence = get_sequence_from_pdb(pdb_file)
    if sequence:
        header = f">sequence_from_{os.path.basename(pdb_file)}"
        fasta_output = f"{header}\n{sequence}\n"
        if output_file:
            with open(output_file, 'w') as f: f.write(fasta_output)
            print(f"Sequência FASTA salva em '{output_file}'")
        else:
            print(fasta_output)

if __name__ == "__main__":
    main()
