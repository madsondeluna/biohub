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
import csv # Importação do módulo para lidar com CSV

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

# --- Banner e Funções Auxiliares ---
def print_banner():
    """Exibe a arte ASCII e informações da ferramenta."""
    # (Implementação do banner mantida)
    banner = r"""
 _._     _,-'""`-._
(,-.`._,'(       |\`-/|
    `-.-' \ )-`( , o o)
          `-    \`_`"'-
█████████████████████████████████████████████████████████████████████
█▌                                                                ▐█
█▌                                                                ▐█
█▌   .______    __    ______    __    __   __   __ .______         ▐█
█▌   |   _  \  |  |  /      \  |  |  |  | |  |  |  | |   _  \        ▐█
█▌   |  |_)  | |  | |  ,----'  |  |__|  | |  |  |  | |  |_)  |       ▐█
█▌   |   _  <  |  | |  |       |   __   | |  |  |  | |   _  <        ▐█
█▌   |  |_)  | |  | |  `----.  |  |  |  | |  `--'  | |  |_)  |       ▐█
█▌   |______/  |__|  \______/  |__|  |__|  \______/  |______/        ▐█
█▌                                                                ▐█
█▌                                                                ▐█
█████████████████████████████████████████████████████████████████████
    """
    print(banner)
    print("Uma Plataforma para Análise de Sequências e Estruturas de Proteínas")
    print("Versão: 1.2.0 | UFMG - Bioinformática")
    print("-" * 65)

def parse_pdb_atoms(pdb_filepath: str):
    # (Implementação mantida)
    atoms = []
    try:
        with open(pdb_filepath, 'r') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom = {
                        "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
                        "res_name": line[17:20].strip(), "res_num": int(line[22:26]),
                        "chain_id": line[21],
                        "element": line[76:78].strip().upper() or line[12:14].strip().upper()
                    }
                    atoms.append(atom)
    except FileNotFoundError:
        print(f"Erro: Arquivo não encontrado em '{pdb_filepath}'", file=sys.stderr)
        return []
    return atoms

def generate_sphere_points(n_points: int):
    # (Implementação mantida)
    points = []
    phi = (1 + math.sqrt(5)) / 2
    for i in range(n_points):
        y = 1 - (2 * i / (n_points - 1))
        radius = math.sqrt(1 - y * y)
        theta = 2 * math.pi * i / phi
        x = math.cos(theta) * radius
        z = math.sin(theta) * radius
        points.append((x, y, z))
    return points

# --- Funções Principais de Análise ---
# (As funções de análise existentes foram mantidas)
def get_sequence_from_pdb(pdb_filepath: str) -> str:
    # (Implementação mantida)
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
        return ""
    return sequence

# (calculate_physicochemical_properties, calculate_intramolecular_contacts, etc., foram mantidas e omitidas para brevidade)
def calculate_physicochemical_properties(sequence: str):
    # (Implementação mantida)
    pass
def calculate_intramolecular_contacts(pdb_filepath: str, threshold: float = 8.0):
    # (Implementação mantida)
    pass
def predict_solvent_exposure(pdb_filepath: str, window_size: int = 9):
    # (Implementação mantida)
    pass
def calculate_sasa(pdb_filepath: str, probe_radius: float, n_points: int):
    # (Implementação mantida)
    pass
def run_apbs_analysis(pdb_filepath: str, no_cleanup: bool):
    # (Implementação mantida)
    pass

# --- NOVAS FUNÇÕES DE CONVERSÃO ---

def handle_pdb_to_fasta(args):
    """Lida com a conversão de PDB para FASTA."""
    sequence = get_sequence_from_pdb(args.pdb_file)
    if sequence:
        header = f">sequence_from_{os.path.basename(args.pdb_file)}"
        fasta_output = f"{header}\n{sequence}\n"
        if args.output:
            with open(args.output, 'w') as f:
                f.write(fasta_output)
            print(f"Sequência FASTA salva em '{args.output}'")
        else:
            print(fasta_output)

def handle_csv_to_fasta(args):
    """Lida com a conversão de CSV para FASTA."""
    try:
        with open(args.csv_file, mode='r', newline='', encoding='utf-8') as infile:
            reader = csv.reader(infile, delimiter=args.delimiter)
            
            id_col_idx, seq_col_idx = -1, -1
            
            if args.header:
                header_row = next(reader)
                try:
                    # Tenta converter para int primeiro (se o usuário passar um número)
                    id_col_idx = int(args.id_col)
                    seq_col_idx = int(args.seq_col)
                except ValueError:
                    # Se falhar, assume que são nomes de colunas
                    if args.id_col in header_row:
                        id_col_idx = header_row.index(args.id_col)
                    if args.seq_col in header_row:
                        seq_col_idx = header_row.index(args.seq_col)
                
                if id_col_idx == -1 or seq_col_idx == -1:
                    print(f"Erro: Coluna de ID ('{args.id_col}') ou Sequência ('{args.seq_col}') não encontrada no cabeçalho.", file=sys.stderr)
                    return
            else:
                try:
                    id_col_idx = int(args.id_col)
                    seq_col_idx = int(args.seq_col)
                except ValueError:
                    print("Erro: Sem cabeçalho, --id-col e --seq-col devem ser índices numéricos.", file=sys.stderr)
                    return

            fasta_records = []
            for i, row in enumerate(reader, 1):
                try:
                    identifier = row[id_col_idx].strip()
                    sequence = row[seq_col_idx].strip().replace(" ", "")
                    if identifier and sequence:
                        fasta_records.append(f">{identifier}\n{sequence}")
                except IndexError:
                    print(f"Aviso: Linha {i+1} no CSV é muito curta e foi ignorada.", file=sys.stderr)
            
            output_content = "\n".join(fasta_records)
            if args.output:
                with open(args.output, 'w') as outfile:
                    outfile.write(output_content)
                print(f"Arquivo FASTA gerado e salvo em '{args.output}'")
            else:
                print(output_content)

    except FileNotFoundError:
        print(f"Erro: Arquivo CSV não encontrado em '{args.csv_file}'", file=sys.stderr)
    except Exception as e:
        print(f"Ocorreu um erro inesperado: {e}", file=sys.stderr)

# --- Configuração da Interface de Linha de Comando ---
def main():
    if len(sys.argv) == 1 or "-h" in sys.argv or "--help" in sys.argv:
        print_banner()

    parser = argparse.ArgumentParser(
        description="BioHub: Uma ferramenta CLI para análise de proteínas.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Função a ser executada")

    # Comando 'fasta' (agora para PDB -> FASTA)
    parser_fasta = subparsers.add_parser("fasta", help="Converte um arquivo PDB em uma sequência FASTA.")
    parser_fasta.add_argument("pdb_file", help="Caminho para o arquivo PDB de entrada.")
    parser_fasta.add_argument("-o", "--output", help="Arquivo de saída (padrão: stdout).")

    # NOVO Comando 'csv2fasta'
    parser_csv = subparsers.add_parser("csv2fasta", help="Converte um arquivo CSV em um formato FASTA.")
    parser_csv.add_argument("csv_file", help="Caminho para o arquivo CSV de entrada.")
    parser_csv.add_argument("-o", "--output", help="Arquivo de saída FASTA (padrão: stdout).")
    parser_csv.add_argument("--id-col", default="0", help="Coluna do identificador (índice ou nome). Padrão: 0.")
    parser_csv.add_argument("--seq-col", default="1", help="Coluna da sequência (índice ou nome). Padrão: 1.")
    parser_csv.add_argument("--header", action="store_true", help="Indica que o CSV tem uma linha de cabeçalho.")
    parser_csv.add_argument("--delimiter", default=",", help="Delimitador do CSV. Padrão: ',' (vírgula).")
    
    # (Restante dos parsers mantidos)
    parser_physchem = subparsers.add_parser("physchem", help="Calcula propriedades físico-químicas de uma sequência.")
    # ... argumentos para physchem
    parser_contacts = subparsers.add_parser("contacts", help="Calcula contatos intramoleculares.")
    # ... argumentos para contacts
    parser_exposure = subparsers.add_parser("exposure", help="Prevê a exposição ao solvente.")
    # ... argumentos para exposure
    parser_sasa = subparsers.add_parser("sasa", help="Calcula a área de superfície acessível ao solvente (SASA).")
    # ... argumentos para sasa
    parser_apbs = subparsers.add_parser("apbs", help="Calcula a energia de solvatação via APBS.")
    # ... argumentos para apbs
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    command_functions = {
        "fasta": handle_pdb_to_fasta,
        "csv2fasta": handle_csv_to_fasta,
        "physchem": lambda a: calculate_physicochemical_properties(a.sequence.upper()),
        "contacts": lambda a: calculate_intramolecular_contacts(a.pdb_file, a.threshold),
        "exposure": lambda a: predict_solvent_exposure(a.pdb_file, a.window),
        "sasa": lambda a: calculate_sasa(a.pdb_file, a.probe_radius, a.num_points),
        "apbs": lambda a: run_apbs_analysis(a.pdb_file, a.no_cleanup)
    }
    command_functions[args.command](args)

if __name__ == "__main__":
    main()
