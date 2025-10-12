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

# --- Funções de Interface ---
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
    print("Versão: 1.6.1 | UFMG - Bioinformática", file=sys.stderr)
    print("-" * 65, file=sys.stderr)

def print_exit_message():
    """Exibe uma mensagem de agradecimento ao final da execução."""
    # Adicionado \n para criar um espaço antes da mensagem
    print("\n" + "=" * 65, file=sys.stderr)
    print("Obrigado por usar o BioHub!", file=sys.stderr)
    print("=" * 65, file=sys.stderr)

# --- Funções Auxiliares (mantidas) ---
def parse_pdb_atoms(pdb_filepath: str):
    # (Implementação mantida)
    pass
def get_sequence_from_pdb(pdb_filepath: str) -> str:
    # (Implementação mantida)
    pass
def write_csv(filepath, header, data_rows):
    # (Implementação mantida)
    pass

# --- Funções de Download, Conversão e Análise (mantidas) ---
def handle_fetch_pdb(args):
    # (Implementação mantida)
    pass
def handle_pdb_to_fasta(args):
    # (Implementação mantida)
    pass
def handle_csv_to_fasta(args):
    # (Implementação mantida)
    pass
def calculate_physicochemical_properties(args):
    # (Implementação mantida)
    pass
def calculate_intramolecular_contacts(args):
    # (Implementação mantida)
    pass
def predict_solvent_exposure(args):
    # (Implementação mantida)
    pass
def calculate_sasa(args):
    # (Implementação mantida)
    pass
def run_apbs_analysis(args):
    # (Implementação mantida)
    pass

# --- Configuração da Interface de Linha de Comando ---
def main():
    if sys.stdout.isatty() or len(sys.argv) == 1:
        print_banner()

    parser = argparse.ArgumentParser(
        description="BioHub: Uma ferramenta CLI para análise de proteínas.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Use: biohub.py [COMANDO] -h para ajuda sobre um comando."
    )
    subparsers = parser.add_subparsers(dest="command", help="Função a ser executada")
    
    # Comandos (configuração mantida)
    parser_fetch = subparsers.add_parser("fetchpdb", help="Baixa um arquivo PDB do RCSB.")
    parser_fetch.add_argument("pdb_id", help="O ID de 4 caracteres do PDB.")
    parser_fetch.add_argument("-o", "--output", help="Nome do arquivo de saída.")

    parser_fasta = subparsers.add_parser("fasta", help="Converte PDB para FASTA.")
    parser_fasta.add_argument("pdb_file", help="Arquivo PDB de entrada.")
    parser_fasta.add_argument("-o", "--output", help="Salva a saída em um arquivo FASTA.")
    
    parser_csv = subparsers.add_parser("csv2fasta", help="Converte CSV para FASTA.")
    parser_csv.add_argument("csv_file", help="Arquivo CSV de entrada.")
    parser_csv.add_argument("-o", "--output", help="Salva a saída em um arquivo FASTA.")
    parser_csv.add_argument("--id-col", default="0", help="Coluna do ID.")
    parser_csv.add_argument("--seq-col", default="1", help="Coluna da sequência.")
    parser_csv.add_argument("--header", action="store_true", help="O CSV contém cabeçalho.")
    parser_csv.add_argument("--delimiter", default=",", help="Delimitador do CSV.")

    analysis_parsers = []
    parser_physchem = subparsers.add_parser("physchem", help="Calcula propriedades físico-químicas.")
    parser_physchem.add_argument("sequence", help="Sequência de aminoácidos.")
    analysis_parsers.append(parser_physchem)
    # (outros parsers de análise omitidos para brevidade)
    
    for p in analysis_parsers:
        p.add_argument("-o", "--output", help="Salva os resultados em um arquivo CSV.")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()
    
    command_functions = {
        "fetchpdb": handle_fetch_pdb,
        "fasta": handle_pdb_to_fasta,
        "csv2fasta": handle_csv_to_fasta,
        "physchem": calculate_physicochemical_properties,
        "contacts": calculate_intramolecular_contacts,
        "exposure": predict_solvent_exposure,
        "sasa": calculate_sasa,
        "apbs": run_apbs_analysis
    }
    
    if args.command in command_functions:
        command_functions[args.command](args)
    else:
        parser.print_help(sys.stderr)

    if sys.stdout.isatty():
        print_exit_message()

if __name__ == "__main__":
    main()
