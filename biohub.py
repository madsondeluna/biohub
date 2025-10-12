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
import shutil
import csv
import urllib.request
from collections import defaultdict

# --- Constantes e Dicionários de Dados ---
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 
    'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 
    'TYR': 'Y', 'VAL': 'V'
}
MOLECULAR_WEIGHT = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16, 'E': 147.13,
    'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17, 'L': 131.17, 'K': 146.19,
    'M': 149.21, 'F': 165.19, 'P': 115.13, 'S': 105.09, 'T': 119.12, 'W': 204.23,
    'Y': 181.19, 'V': 117.15
}
KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5, 'Q': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}
PKA_VALUES = {
    'C-term': 3.65, 'N-term': 8.0, 'D': 3.9, 'E': 4.07, 'H': 6.5, 'C': 8.5,
    'Y': 10.0, 'K': 10.0, 'R': 12.0
}
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80, 'F': 1.47,
    'CL': 1.75, 'BR': 1.85, 'I': 1.98, 'DEFAULT': 1.70
}
DIWV = {
    'A': {'A': 42,'R':-24,'N':-36,'D':-52,'C':  8,'Q':-38,'E':-43,'G': 60,'H':-24,'I': 27,'L': 15,'K':-35,'M':-15,'F':-16,'P':-13,'S': 36,'T': 21,'W':-33,'Y':-15,'V': 42},
    'R': {'A':-24,'R': 42,'N':-15,'D':-60,'C':-44,'Q': -3,'E':-27,'G':-41,'H':  9,'I':-29,'L':-33,'K': 16,'M': -9,'F':-44,'P':  3,'S':-12,'T':-21,'W': 27,'Y':-33,'V':-36},
    'N': {'A':-36,'R':-15,'N': 42,'D': 21,'C':-35,'Q': 13,'E': -6,'G':  6,'H': -2,'I':-33,'L':-44,'K':  1,'M':-42,'F':-33,'P':-42,'S':  4,'T': -8,'W':-29,'Y':-15,'V':-44},
    'D': {'A':-52,'R':-60,'N': 21,'D': 42,'C':-60,'Q': -6,'E': 36,'G':-13,'H':-24,'I':-60,'L':-60,'K':-21,'M':-60,'F':-60,'P':-33,'S': -6,'T':-27,'W':-60,'Y':-42,'V':-60},
    'C': {'A':  8,'R':-44,'N':-35,'D':-60,'C': 42,'Q':-42,'E':-60,'G':-29,'H':-29,'I': 13,'L': -6,'K':-60,'M': -8,'F':-13,'P':-41,'S': -3,'T': -3,'W':-44,'Y': 24,'V': 12},
    'Q': {'A':-38,'R': -3,'N': 13,'D': -6,'C':-42,'Q': 42,'E': 16,'G':-13,'H': 21,'I':-33,'L':-27,'K':  4,'M': -6,'F':-42,'P': -9,'S':-13,'T':-16,'W':-21,'Y':-27,'V':-38},
    'E': {'A':-43,'R':-27,'N': -6,'D': 36,'C':-60,'Q': 16,'E': 42,'G':-21,'H': -8,'I':-42,'L':-41,'K': -1,'M':-29,'F':-42,'P':-24,'S': -1,'T':-13,'W':-41,'Y':-27,'V':-42},
    'G': {'A': 60,'R':-41,'N':  6,'D':-13,'C':-29,'Q':-13,'E':-21,'G': 42,'H':-44,'I':-44,'L':-44,'K':-24,'M':-44,'F':-44,'P':-21,'S': 16,'T': -6,'W':-44,'Y':-44,'V':-44},
    'H': {'A':-24,'R':  9,'N': -2,'D':-24,'C':-29,'Q': 21,'E': -8,'G':-44,'H': 42,'I':-24,'L':-29,'K': -8,'M':-15,'F':-16,'P': -8,'S':-16,'T':-16,'W': -3,'Y':  0,'V':-29},
    'I': {'A': 27,'R':-29,'N':-33,'D':-60,'C': 13,'Q':-33,'E':-42,'G':-44,'H':-24,'I': 42,'L': 36,'K':-33,'M': 24,'F': 13,'P':-33,'S':-21,'T':  8,'W':-27,'Y': -1,'V': 36},
    'L': {'A': 15,'R':-33,'N':-44,'D':-60,'C': -6,'Q':-27,'E':-41,'G':-44,'H':-29,'I': 36,'L': 42,'K':-33,'M': 33,'F': 21,'P':-33,'S':-33,'T': -9,'W':-13,'Y': -8,'V': 33},
    'K': {'A':-35,'R': 16,'N':  1,'D':-21,'C':-60,'Q':  4,'E': -1,'G':-24,'H': -8,'I':-33,'L':-33,'K': 42,'M':  0,'F':-42,'P':-16,'S': -1,'T': -7,'W':-21,'Y':-21,'V':-38},
    'M': {'A':-15,'R': -9,'N':-42,'D':-60,'C': -8,'Q': -6,'E':-29,'G':-44,'H':-15,'I': 24,'L': 33,'K':  0,'M': 42,'F': 13,'P':-29,'S':-24,'T': -9,'W':-21,'Y': -8,'V': 21},
    'F': {'A':-16,'R':-44,'N':-33,'D':-60,'C':-13,'Q':-42,'E':-42,'G':-44,'H':-16,'I': 13,'L': 21,'K':-42,'M': 13,'F': 42,'P':-44,'S':-27,'T':-21,'W': 21,'Y': 36,'V': -1},
    'P': {'A':-13,'R':  3,'N':-42,'D':-33,'C':-41,'Q': -9,'E':-24,'G':-21,'H': -8,'I':-33,'L':-33,'K':-16,'M':-29,'F':-44,'P': 42,'S':  1,'T': -9,'W':-44,'Y':-33,'V':-27},
    'S': {'A': 36,'R':-12,'N':  4,'D': -6,'C': -3,'Q':-13,'E': -1,'G': 16,'H':-16,'I':-21,'L':-33,'K': -1,'M':-24,'F':-27,'P':  1,'S': 42,'T': 16,'W':-24,'Y':-15,'V':-21},
    'T': {'A': 21,'R':-21,'N': -8,'D':-27,'C': -3,'Q':-16,'E':-13,'G': -6,'H':-16,'I':  8,'L': -9,'K': -7,'M': -9,'F':-21,'P': -9,'S': 16,'T': 42,'W':-29,'Y':-16,'V':  4},
    'W': {'A':-33,'R': 27,'N':-29,'D':-60,'C':-44,'Q':-21,'E':-41,'G':-44,'H': -3,'I':-27,'L':-13,'K':-21,'M':-21,'F': 21,'P':-44,'S':-24,'T':-29,'W': 42,'Y': 21,'V':-29},
    'Y': {'A':-15,'R':-33,'N':-15,'D':-42,'C': 24,'Q':-27,'E':-27,'G':-44,'H':  0,'I': -1,'L': -8,'K':-21,'M': -8,'F': 36,'P':-33,'S':-15,'T':-16,'W': 21,'Y': 42,'V':-13},
    'V': {'A': 42,'R':-36,'N':-44,'D':-60,'C': 12,'Q':-38,'E':-42,'G':-44,'H':-29,'I': 36,'L': 33,'K':-38,'M': 21,'F': -1,'P':-27,'S':-21,'T':  4,'W':-29,'Y':-13,'V': 42}
}

# --- Funções de Interface e Auxiliares ---
def print_banner():
    """Exibe a arte ASCII e informações da ferramenta."""
    banner = r"""

 _._     _,-'""`-._
(,-.`._,'(       |\`-/|
    `-.-' \ )-`( , o o)
          `-    \`_`"'-
█████████████████████████████████████████████████████████████████████
█▌                                                                 ▐█
█▌                                                                 ▐█
█▌    .______    __    ______    __    __   __    __  .______      ▐█
█▌    |   _  \  |  |  /  __  \  |  |  |  | |  |  |  | |   _  \     ▐█
█▌    |  |_)  | |  | |  |  |  | |  |__|  | |  |  |  | |  |_)  |    ▐█
█▌    |   _  <  |  | |  |  |  | |   __   | |  |  |  | |   _  <     ▐█
█▌    |  |_)  | |  | |  `--'  | |  |  |  | |  `--'  | |  |_)  |    ▐█
█▌    |______/  |__|  \______/  |__|  |__|  \______/  |______/     ▐█
█▌                                                                 ▐█
█▌                                                                 ▐█
█████████████████████████████████████████████████████████████████████

    """
    print(banner, file=sys.stderr)
    print("Uma Plataforma para Análise de Sequências e Estruturas de Proteínas", file=sys.stderr)
    print("Versão: 1.9.0 | UFMG - Bioinformática", file=sys.stderr)
    print("Autores: ACDS, AKNNA, LSRS, LHDS, MADLA", file=sys.stderr)
    print("-" * 65, file=sys.stderr)

def print_exit_message():
    """Exibe uma mensagem de agradecimento ao final da execução."""
    print("\n" + "=" * 65, file=sys.stderr)
    print("Obrigado por usar o BioHub! Essa aplicação foi feita com <3 e cafê..." + "\nAh, se vocè ver um gatinho carente no ICB, faça carinho nele!", file=sys.stderr)
    print("=" * 65, file=sys.stderr)

def parse_pdb_atoms(pdb_filepath: str):
    atoms = []
    try:
        with open(pdb_filepath, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atoms.append({
                        "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
                        "res_name": line[17:20].strip(), "res_num": int(line[22:26]),
                        "chain_id": line[21],
                        "element": line[76:78].strip().upper() or line[12:14].strip().upper()
                    })
    except FileNotFoundError:
        print(f"Erro: Arquivo não encontrado em '{pdb_filepath}'", file=sys.stderr)
    return atoms

def get_sequence_from_pdb(pdb_filepath: str) -> str:
    sequence = ""
    processed_residues = set()
    first_chain_id = None
    try:
        with open(pdb_filepath, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    current_chain_id = line[21]
                    if first_chain_id is None: first_chain_id = current_chain_id
                    if current_chain_id == first_chain_id:
                        res_id = (int(line[22:26]), line[21])
                        if res_id not in processed_residues:
                            res_name = line[17:20]
                            if res_name in THREE_TO_ONE:
                                sequence += THREE_TO_ONE[res_name]
                                processed_residues.add(res_id)
    except FileNotFoundError:
        print(f"Erro: Arquivo não encontrado em '{pdb_filepath}'", file=sys.stderr)
    return sequence

def extract_pdb_header_info(pdb_filepath: str):
    """Extrai informações do cabeçalho de um arquivo PDB."""
    info = defaultdict(str)
    try:
        with open(pdb_filepath, 'r') as f:
            for line in f:
                if line.startswith("HEADER"):
                    info['classification'] = line[10:50].strip()
                    info['dep_date'] = line[50:59].strip()
                elif line.startswith("TITLE"):
                    info['title'] += line[10:80].strip() + " "
                elif line.startswith("SOURCE"):
                    if "ORGANISM_SCIENTIFIC" in line:
                        info['organism'] += line.split(":")[-1].strip().replace(';','') + " "
    except Exception as e:
        print(f"Aviso: Não foi possível extrair informações do cabeçalho do PDB: {e}", file=sys.stderr)
    
    # Limpa espaços extras
    for key in info:
        info[key] = ' '.join(info[key].split())
    return info
    
def write_csv(filepath, header, data_rows):
    try:
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(data_rows)
        print(f"Resultados salvos com sucesso em '{filepath}'", file=sys.stderr)
    except IOError as e:
        print(f"Erro ao salvar o arquivo: {e}", file=sys.stderr)

def generate_sphere_points(n_points: int):
    points = []
    if n_points <= 0: return points
    phi = (1 + math.sqrt(5)) / 2
    for i in range(n_points):
        y = 1 - (2 * i / (n_points - 1)) if n_points > 1 else 0
        radius = math.sqrt(1 - y * y)
        theta = 2 * math.pi * i / phi
        points.append((math.cos(theta) * radius, y, math.sin(theta) * radius))
    return points

# --- Funções de Download, Conversão e Análise ---
def handle_fetch_pdb(args):
    pdb_id = args.pdb_id.upper()
    if len(pdb_id) != 4:
        print(f"Erro: O ID do PDB '{pdb_id}' é inválido. Deve conter 4 caracteres.", file=sys.stderr)
        return
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_file = args.output if args.output else f"{pdb_id}.pdb"
    try:
        print(f"Baixando {pdb_id} de {url}...", file=sys.stderr)
        with urllib.request.urlopen(url) as response, open(output_file, 'wb') as out_file:
            data = response.read()
            if not data: raise ValueError("Arquivo PDB recebido está vazio.")
            out_file.write(data)
        print(f"Arquivo PDB salvo com sucesso em '{output_file}'", file=sys.stderr)
        
        # Exibe informações do cabeçalho
        info = extract_pdb_header_info(output_file)
        if info:
            print("\n--- Informações da Estrutura ---", file=sys.stderr)
            if info['title']: print(f"  Título    : {info['title']}", file=sys.stderr)
            if info['dep_date']: print(f"  Data      : {info['dep_date']}", file=sys.stderr)
            if info['classification']: print(f"  Classe    : {info['classification']}", file=sys.stderr)
            if info['organism']: print(f"  Organismo : {info['organism']}", file=sys.stderr)

    except urllib.error.HTTPError as e:
        print(f"Erro ao baixar o arquivo: Não foi possível encontrar o PDB ID '{pdb_id}'. (HTTP {e.code})", file=sys.stderr)
        if os.path.exists(output_file): os.remove(output_file)
    except Exception as e:
        print(f"Ocorreu um erro inesperado: {e}", file=sys.stderr)
        if os.path.exists(output_file): os.remove(output_file)

def handle_pdb_to_fasta(args):
    sequence = get_sequence_from_pdb(args.pdb_file)
    if sequence:
        header = f">sequence_from_{os.path.basename(args.pdb_file)}"
        fasta_output = f"{header}\n{sequence}\n"
        if args.output:
            with open(args.output, 'w') as f: f.write(fasta_output)
            print(f"Sequência FASTA salva em '{args.output}'", file=sys.stderr)
        else:
            print(fasta_output)

def handle_csv_to_fasta(args):
    try:
        with open(args.csv_file, mode='r', newline='', encoding='utf-8') as infile:
            reader = csv.reader(infile, delimiter=args.delimiter)
            id_col_idx, seq_col_idx = -1, -1
            if args.header:
                header_row = next(reader)
                try:
                    id_col_idx, seq_col_idx = int(args.id_col), int(args.seq_col)
                except ValueError:
                    if args.id_col in header_row: id_col_idx = header_row.index(args.id_col)
                    if args.seq_col in header_row: seq_col_idx = header_row.index(args.seq_col)
                if id_col_idx == -1 or seq_col_idx == -1:
                    print(f"Erro: Coluna de ID ('{args.id_col}') ou Sequência ('{args.seq_col}') não encontrada.", file=sys.stderr)
                    return
            else:
                try:
                    id_col_idx, seq_col_idx = int(args.id_col), int(args.seq_col)
                except ValueError:
                    print("Erro: Sem cabeçalho, as colunas devem ser índices numéricos.", file=sys.stderr)
                    return
            fasta_records = [f">{row[id_col_idx].strip()}\n{row[seq_col_idx].strip().replace(' ', '')}" for row in reader if row and len(row) > max(id_col_idx, seq_col_idx)]
            output_content = "\n".join(fasta_records)
            if args.output:
                with open(args.output, 'w') as outfile: outfile.write(output_content)
                print(f"Arquivo FASTA salvo em '{args.output}'", file=sys.stderr)
            else:
                print(output_content)
    except FileNotFoundError:
        print(f"Erro: Arquivo CSV não encontrado em '{args.csv_file}'", file=sys.stderr)
    except Exception as e:
        print(f"Ocorreu um erro: {e}", file=sys.stderr)

def calculate_physicochemical_properties(args):
    sequence = args.sequence.upper()
    if not sequence:
        print("Erro: Sequência de entrada está vazia.", file=sys.stderr)
        return
    length = len(sequence)
    aa_composition = {aa: sequence.count(aa) for aa in MOLECULAR_WEIGHT.keys()}
    mw = sum(MOLECULAR_WEIGHT.get(aa, 0) * count for aa, count in aa_composition.items()) - (length - 1) * 18.015
    gravy = sum(KYTE_DOOLITTLE.get(aa, 0) for aa in sequence) / length if length > 0 else 0
    pi = 0.0
    min_charge_diff = float('inf')
    for ph_int in range(1401):
        ph = ph_int * 0.01
        net_charge = (10**PKA_VALUES['N-term']) / (10**PKA_VALUES['N-term'] + 10**ph) - (10**ph) / (10**PKA_VALUES['C-term'] + 10**ph)
        for aa, count in aa_composition.items():
            if aa in ['R', 'H', 'K']: net_charge += count * (10**PKA_VALUES[aa]) / (10**PKA_VALUES[aa] + 10**ph)
            elif aa in ['D', 'E', 'C', 'Y']: net_charge -= count * (10**ph) / (10**PKA_VALUES[aa] + 10**ph)
        if abs(net_charge) < min_charge_diff:
            min_charge_diff, pi = abs(net_charge), ph
    a, b = 2.9, 3.9
    aliphatic_index = (aa_composition.get('A', 0) + a * aa_composition.get('V', 0) + b * (aa_composition.get('I', 0) + aa_composition.get('L', 0))) / length * 100 if length > 0 else 0
    instability_index = (10 / length) * sum(DIWV[sequence[i]][sequence[i+1]] for i in range(length - 1)) if length > 1 else 0
    stability = "Estável" if instability_index < 40 else "Instável"
    n_term = sequence[0]
    half_life_mammal = ">10 horas" if n_term in ['A','C','G','M','P','S','T','V'] else "2-30 min" if n_term in ['I','L','F','W','Y','D','E','N','Q'] else "Desconhecido"
    def count_group(group_aas): return sum(aa_composition.get(aa, 0) for aa in group_aas)
    acidic, basic = count_group({'D', 'E'}), count_group({'R', 'K', 'H'})
    polar, non_polar = count_group({'N', 'Q', 'S', 'T', 'Y', 'C'}), count_group({'A', 'V', 'L', 'I', 'P', 'F', 'W', 'M', 'G'})
    results = [
        ("Comprimento", length), ("Peso Molecular (Da)", f"{mw:.2f}"),
        ("Ponto Isoelétrico (pI)", f"{pi:.2f}"), ("GRAVY (Hidropaticidade)", f"{gravy:.3f}"),
        ("Índice Alifático", f"{aliphatic_index:.2f}"), ("Índice de Instabilidade", f"{instability_index:.2f} ({stability})"),
        ("Meia-Vida (Mamíferos, in vitro)", half_life_mammal), ("Total de Resíduos Ácidos (Asp+Glu)", acidic),
        ("Total de Resíduos Básicos (Arg+Lys+His)", basic), ("Total de Resíduos Polares", polar),
        ("Total de Resíduos Apolares", non_polar)
    ]
    if args.output:
        write_csv(args.output, ["Propriedade", "Valor"], results)
    else:
        print("--- Propriedades Físico-Químicas ---")
        for key, value in results: print(f"{key:<35}: {value}")
        print("\n--- Composição de Aminoácidos (%) ---")
        for aa, count in sorted(aa_composition.items()):
            if count > 0: print(f"{aa}: {count:<5} ({(count/length)*100:.2f}%)", end="  ")
        print()

def calculate_intramolecular_contacts(args):
    ca_atoms = {}
    if not os.path.exists(args.pdb_file):
        print(f"Erro: Arquivo não encontrado em '{args.pdb_file}'", file=sys.stderr)
        return
    with open(args.pdb_file, 'r') as f:
        first_chain_id = None
        for line in f:
            if line.startswith("ATOM") and line[13:16].strip() == "CA":
                if first_chain_id is None: first_chain_id = line[21]
                if line[21] == first_chain_id:
                    ca_atoms[int(line[22:26])] = tuple(float(line[i:i+8]) for i in [30, 38, 46])
    residues = sorted(ca_atoms.keys())
    contacts = [(r1, r2, round(math.sqrt(sum((c1-c2)**2 for c1,c2 in zip(ca_atoms[r1], ca_atoms[r2]))), 3))
                for i, r1 in enumerate(residues) for r2 in residues[i+1:]
                if abs(r1 - r2) > 1 and math.sqrt(sum((c1-c2)**2 for c1,c2 in zip(ca_atoms[r1], ca_atoms[r2]))) <= args.threshold]
    if args.output:
        write_csv(args.output, ["Residuo_1", "Residuo_2", "Distancia_A"], contacts)
    else:
        print(f"--- Contatos Intramoleculares (Limiar = {args.threshold:.1f} Å) ---")
        if not contacts: print("Nenhum contato encontrado.")
        else:
            for c in contacts: print(f"Res {c[0]} - Res {c[1]}: {c[2]:.3f} Å")

def predict_solvent_exposure(args):
    sequence = get_sequence_from_pdb(args.pdb_file)
    if not sequence: return
    results = []
    half_window = args.window // 2
    for i in range(len(sequence)):
        window_seq = sequence[max(0, i - half_window) : min(len(sequence), i + half_window + 1)]
        score = sum(KYTE_DOOLITTLE.get(aa, 0) for aa in window_seq) / len(window_seq)
        results.append([i + 1, sequence[i], f"{score:.3f}"])
    if args.output:
        write_csv(args.output, ["Posicao", "Residuo", "Score_Hidropatia"], results)
    else:
        print(f"--- Predição de Exposição (Kyte-Doolittle, Janela={args.window}) ---")
        print("Pos. | Res. | Score de Hidropatia")
        for row in results: print(f"{row[0]:<4} | {row[1]:<4} | {row[2]}")

def calculate_sasa(args):
    atoms = parse_pdb_atoms(args.pdb_file)
    if not atoms: return
    sphere_points = generate_sphere_points(args.num_points)
    total_sasa = 0.0
    sasa_per_residue = {}
    for i, atom_i in enumerate(atoms):
        radius_i = VDW_RADII.get(atom_i["element"], VDW_RADII['DEFAULT'])
        extended_radius = radius_i + args.probe_radius
        accessible_points = 0
        for sp in sphere_points:
            point_is_accessible = True
            point = (atom_i["x"] + extended_radius * sp[0], atom_i["y"] + extended_radius * sp[1], atom_i["z"] + extended_radius * sp[2])
            for j, atom_j in enumerate(atoms):
                if i == j: continue
                radius_j = VDW_RADII.get(atom_j["element"], VDW_RADII['DEFAULT'])
                if sum((p-c)**2 for p,c in zip(point, (atom_j["x"], atom_j["y"], atom_j["z"]))) < (radius_j + args.probe_radius)**2:
                    point_is_accessible = False; break
            if point_is_accessible: accessible_points += 1
        atom_sasa = (accessible_points / args.num_points) * 4.0 * math.pi * extended_radius**2 if args.num_points > 0 else 0
        total_sasa += atom_sasa
        res_id = (atom_i["res_num"], atom_i["res_name"])
        sasa_per_residue[res_id] = sasa_per_residue.get(res_id, 0) + atom_sasa
    print(f"SASA Total da Molécula: {total_sasa:.2f} Å²", file=sys.stderr)
    results_data = [[f"{res_id[1]} {res_id[0]}", f"{sasa:.2f}"] for res_id, sasa in sorted(sasa_per_residue.items())]
    if args.output:
        write_csv(args.output, ["Residuo", "SASA_A2"], results_data)
    else:
        print("--- SASA por Resíduo ---")
        print("Resíduo  | SASA (Å²)")
        for row in results_data: print(f"{row[0]:<8} | {row[1]}")

def run_apbs_analysis(args):
    for exe in ["pdb2pqr", "apbs"]:
        if not shutil.which(exe):
            print(f"Erro: '{exe}' não encontrado no PATH.", file=sys.stderr); return
    work_dir = tempfile.mkdtemp()
    try:
        pqr_path = os.path.join(work_dir, f"{os.path.basename(args.pdb_file)}.pqr")
        print(f"1. Executando pdb2pqr...", file=sys.stderr)
        subprocess.run(["pdb2pqr", "--ff=amber", args.pdb_file, pqr_path], check=True, capture_output=True, text=True, timeout=60)
        apbs_in_path = os.path.join(work_dir, "apbs.in")
        with open(apbs_in_path, 'w') as f: f.write(f"read\n    mol pqr {os.path.basename(pqr_path)}\nend\nelec\n    mg-auto\nend\nquit\n")
        print("2. Executando APBS...", file=sys.stderr)
        result = subprocess.run(["apbs", "apbs.in"], check=True, capture_output=True, text=True, cwd=work_dir, timeout=300)
        print("3. Analisando resultados...", file=sys.stderr)
        energy = "Não encontrada"
        for line in result.stdout.splitlines():
            if "Global net ELEC energy" in line:
                energy = f"{float(line.split()[-2]):.4f} kJ/mol"; break
        print(f"--- Energia de Solvatação Eletrostática ---")
        print(f"Energia: {energy}")
    except subprocess.CalledProcessError as e:
        print(f"Erro ao executar processo externo:\n{e.stderr}", file=sys.stderr)
    finally:
        if not args.no_cleanup: shutil.rmtree(work_dir)
        else: print(f"Arquivos intermediários mantidos em: {work_dir}", file=sys.stderr)

# --- Configuração da Interface de Linha de Comando ---
def main():
    if sys.stdout.isatty() or len(sys.argv) == 1:
        print_banner()

    parser = argparse.ArgumentParser(
        description="BioHub: Uma ferramenta CLI para análise de proteínas.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Use: biohub.py [COMANDO] -h para ajuda detalhada sobre um comando."
    )
    subparsers = parser.add_subparsers(dest="command", help="\nComando a ser executado.")

    
    # --- fetchpdb ---
    parser_fetch = subparsers.add_parser("fetchpdb", help="Baixa um arquivo PDB do RCSB.", formatter_class=argparse.RawTextHelpFormatter)
    parser_fetch.add_argument("pdb_id", metavar="PDB_ID", help="O ID de 4 caracteres do PDB a ser baixado (ex: 1A2B).")
    parser_fetch.add_argument("-o", "--output", metavar="ARQUIVO", help="Nome do arquivo de saída (padrão: [PDB_ID].pdb).")

    # --- fasta ---
    parser_fasta = subparsers.add_parser("fasta", help="Converte um arquivo PDB em uma sequência FASTA.", formatter_class=argparse.RawTextHelpFormatter)
    parser_fasta.add_argument("pdb_file", metavar="ARQUIVO_PDB", help="Caminho para o arquivo PDB de entrada.")
    parser_fasta.add_argument("-o", "--output", metavar="ARQUIVO", help="Salva a saída em um arquivo FASTA (padrão: stdout).")
    
    # --- csv2fasta ---
    parser_csv = subparsers.add_parser("csv2fasta", help="Converte um arquivo CSV em um formato FASTA.", formatter_class=argparse.RawTextHelpFormatter)
    parser_csv.add_argument("csv_file", metavar="ARQUIVO_CSV", help="Caminho para o arquivo CSV de entrada.")
    parser_csv.add_argument("-o", "--output", metavar="ARQUIVO", help="Salva a saída em um arquivo FASTA (padrão: stdout).")
    parser_csv.add_argument("--id-col", metavar="COLUNA", default="0", help="Coluna do identificador (nome ou índice baseado em 0). Padrão: 0.")
    parser_csv.add_argument("--seq-col", metavar="COLUNA", default="1", help="Coluna da sequência (nome ou índice baseado em 0). Padrão: 1.")
    parser_csv.add_argument("--header", action="store_true", help="Flag para indicar que a primeira linha do CSV é um cabeçalho.")
    parser_csv.add_argument("--delimiter", metavar="CHAR", default=",", help="Caractere usado como delimitador no CSV. Padrão: ','.")

    # --- physchem ---
    parser_physchem = subparsers.add_parser("physchem", help="Calcula um conjunto expandido de propriedades físico-químicas.", formatter_class=argparse.RawTextHelpFormatter)
    parser_physchem.add_argument("sequence", metavar="SEQUENCIA", help="A sequência de aminoácidos a ser analisada.")
    parser_physchem.add_argument("-o", "--output", metavar="ARQUIVO_CSV", help="Salva os resultados em um arquivo CSV.")

    # --- contacts ---
    parser_contacts = subparsers.add_parser("contacts", help="Calcula contatos intramoleculares com base na distância entre C-Alfas.", formatter_class=argparse.RawTextHelpFormatter)
    parser_contacts.add_argument("pdb_file", metavar="ARQUIVO_PDB", help="Caminho para o arquivo PDB de entrada.")
    parser_contacts.add_argument("-t", "--threshold", metavar="FLOAT", type=float, default=8.0, help="Distância máxima em Angstroms para considerar um contato. Padrão: 8.0.")
    parser_contacts.add_argument("-o", "--output", metavar="ARQUIVO_CSV", help="Salva os resultados em um arquivo CSV.")

    # --- exposure ---
    parser_exposure = subparsers.add_parser("exposure", help="Prevê regiões de exposição ao solvente com a escala Kyte-Doolittle.", formatter_class=argparse.RawTextHelpFormatter)
    parser_exposure.add_argument("pdb_file", metavar="ARQUIVO_PDB", help="Caminho para o arquivo PDB de entrada.")
    parser_exposure.add_argument("-w", "--window", metavar="INT", type=int, default=9, help="Tamanho da janela deslizante para o cálculo da média. Deve ser ímpar. Padrão: 9.")
    parser_exposure.add_argument("-o", "--output", metavar="ARQUIVO_CSV", help="Salva os resultados em um arquivo CSV.")
    
    # --- sasa ---
    parser_sasa = subparsers.add_parser("sasa", help="Calcula a Área de Superfície Acessível ao Solvente (SASA).", formatter_class=argparse.RawTextHelpFormatter)
    parser_sasa.add_argument("pdb_file", metavar="ARQUIVO_PDB", help="Caminho para o arquivo PDB de entrada.")
    parser_sasa.add_argument("--probe-radius", metavar="FLOAT", type=float, default=1.4, help="Raio da sonda do solvente em Angstroms (padrão: 1.4 para água).")
    parser_sasa.add_argument("--num-points", metavar="INT", type=int, default=960, help="Número de pontos na superfície de cada átomo para o cálculo. Padrão: 960.")
    parser_sasa.add_argument("-o", "--output", metavar="ARQUIVO_CSV", help="Salva os resultados por resíduo em um arquivo CSV.")

    # --- apbs ---
    parser_apbs = subparsers.add_parser("apbs", help="Calcula a energia de solvatação eletrostática (requer PDB2PQR e APBS).", formatter_class=argparse.RawTextHelpFormatter)
    parser_apbs.add_argument("pdb_file", metavar="ARQUIVO_PDB", help="Caminho para o arquivo PDB de entrada.")
    parser_apbs.add_argument("--no-cleanup", action="store_true", help="Previne a remoção dos arquivos temporários (PQR, apbs.in, etc.).")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()
    
    command_functions = {
        "fetchpdb": handle_fetch_pdb, "fasta": handle_pdb_to_fasta, "csv2fasta": handle_csv_to_fasta,
        "physchem": calculate_physicochemical_properties, "contacts": calculate_intramolecular_contacts,
        "exposure": predict_solvent_exposure, "sasa": calculate_sasa, "apbs": run_apbs_analysis
    }
    
    if args.command in command_functions:
        command_functions[args.command](args)
    else:
        parser.print_help(sys.stderr)

    if sys.stdout.isatty():
        print_exit_message()

if __name__ == "__main__":
    main()

