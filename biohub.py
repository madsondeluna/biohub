#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioHub: Uma ferramenta de linha de comando para análises bioinformáticas,
que criei para ter o mínimo de dependências externas.
"""

# Importações de Módulos
# Módulos padrão do Python que vou precisar.
import sys          # Para interagir com o sistema, como ler argumentos e escrever em stderr.
import math         # Para funções matemáticas como raiz quadrada (sqrt) e constantes como pi.
import argparse     # Essencial para criar a interface de linha de comando (CLI) de forma organizada.
import subprocess   # Para executar programas externos, como o PDB2PQR e o APBS.
import os           # Para interagir com o sistema operacional, como manipular nomes de arquivos e caminhos.
import tempfile     # Para criar diretórios temporários para os arquivos do APBS.
import shutil       # Para remover os diretórios temporários.
import csv          # Para ler e escrever arquivos no formato CSV.
import urllib.request # Para baixar arquivos da internet, especificamente os PDBs do RCSB.
from collections import defaultdict # Um tipo de dicionário que cria um valor padrão para chaves que ainda não existem.

# Importação opcional do módulo de visualização
try:
    import biohub_viz
    HAS_VIZ = True
except ImportError:
    HAS_VIZ = False

# Constantes e Dicionários de Dados 
# Dicionário para converter o código de 3 letras de aminoácidos para 1 letra.
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 
    'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 
    'TYR': 'Y', 'VAL': 'V'
}
# Peso molecular de cada aminoácido em Daltons (Da).
MOLECULAR_WEIGHT = {
    'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16, 'E': 147.13,
    'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17, 'L': 131.17, 'K': 146.19,
    'M': 149.21, 'F': 165.19, 'P': 115.13, 'S': 105.09, 'T': 119.12, 'W': 204.23,
    'Y': 181.19, 'V': 117.15
}
# Escala de hidropaticidade de Kyte-Doolittle. Valores mais altos indicam maior hidrofobicidade.
KYTE_DOOLITTLE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5, 'Q': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}
# Valores de pKa para os grupos ionizáveis, usados para calcular o ponto isoelétrico (pI).
PKA_VALUES = {
    'C-term': 3.65, 'N-term': 8.0, 'D': 3.9, 'E': 4.07, 'H': 6.5, 'C': 8.5,
    'Y': 10.0, 'K': 10.0, 'R': 12.0
}
# Raios de Van der Waals (em Angstroms) para diferentes elementos, usados no cálculo do SASA.
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80, 'F': 1.47,
    'CL': 1.75, 'BR': 1.85, 'I': 1.98, 'DEFAULT': 1.70 # Um valor padrão para átomos não listados.
}
# Matriz de pesos (DIWV) usada para calcular o índice de instabilidade de uma proteína.
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

# Funções de Interface e Auxiliares

def print_banner():
    """Esta função simplesmente exibe a arte ASCII e informações da ferramenta no início da execução."""
    banner = r"""
    
    ██████╗ ██╗ ██████╗ ██╗  ██╗██╗   ██╗██████╗ 
    ██╔══██╗██║██╔═══██╗██║  ██║██║   ██║██╔══██╗
    ██████╔╝██║██║   ██║███████║██║   ██║██████╔╝
    ██╔══██╗██║██║   ██║██╔══██║██║   ██║██╔══██╗
    ██████╔╝██║╚██████╔╝██║  ██║╚██████╔╝██████╔╝
    ╚═════╝ ╚═╝ ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ 

    """
    # Uso sys.stderr para que a interface não seja misturada com a saída de dados (stdout).
    print(banner, file=sys.stderr)
    print("Uma Plataforma para Análise de Sequências e Estruturas de Proteínas", file=sys.stderr)
    print("Versão: 0.2.0 | UFMG - Bioinformática", file=sys.stderr)
    print("Autores: ACDS, AKNNA, LSRS, LHS & MADLA", file=sys.stderr)
    print("-" * 68, file=sys.stderr)

def print_exit_message():
    """Exibe uma mensagem de agradecimento ao final da execução..."""
    print("\n" + "=" * 68, file=sys.stderr)
    print("Obrigado por usar o BioHub! Essa aplicação foi feita com <3 e café...", file=sys.stderr)
    print("=" * 68, file=sys.stderr)

def parse_pdb_atoms(pdb_filepath: str):
    """Lê um arquivo PDB e extrai as coordenadas e informações de cada átomo."""
    atoms = []
    try:
        with open(pdb_filepath, 'r') as f:
            for line in f:
                # Linhas ATOM e HETATM contêm as informações dos átomos.
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Extraio as informações de cada coluna, convertendo para o tipo correto.
                    atoms.append({
                        "atom_num": int(line[6:11]),  # Número serial do átomo
                        "atom_name": line[12:16].strip(),  # Nome do átomo (ex: CA, CB, N)
                        "res_name": line[17:20].strip(), 
                        "chain_id": line[21],
                        "res_num": int(line[22:26]),
                        "x": float(line[30:38]), 
                        "y": float(line[38:46]), 
                        "z": float(line[46:54]),
                        # Tento pegar o elemento da coluna 76-78; se não tiver, pego da 12-14.
                        "element": line[76:78].strip().upper() or line[12:14].strip().upper()
                    })
    except FileNotFoundError:
        print(f"Erro: Arquivo não encontrado em '{pdb_filepath}'", file=sys.stderr)
    return atoms

def get_sequence_from_pdb(pdb_filepath: str) -> str:
    """Extrai a sequência de aminoácidos da primeira cadeia de um arquivo PDB."""
    sequence = ""
    processed_residues = set() # Para não adicionar o mesmo resíduo várias vezes.
    first_chain_id = None # Vou pegar a sequência apenas da primeira cadeia que encontrar.
    try:
        with open(pdb_filepath, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    current_chain_id = line[21]
                    # Define a primeira cadeia encontrada como a cadeia de referência.
                    if first_chain_id is None: first_chain_id = current_chain_id
                    
                    # Processo a linha apenas se ela pertencer à primeira cadeia.
                    if current_chain_id == first_chain_id:
                        res_id = (int(line[22:26]), line[21]) # Identificador único do resíduo (número, cadeia).
                        if res_id not in processed_residues:
                            res_name = line[17:20]
                            # Converte o nome de 3 letras para 1 letra e adiciona à sequência.
                            if res_name in THREE_TO_ONE:
                                sequence += THREE_TO_ONE[res_name]
                                processed_residues.add(res_id)
    except FileNotFoundError:
        print(f"Erro: Arquivo não encontrado em '{pdb_filepath}'", file=sys.stderr)
    return sequence

def extract_pdb_header_info(pdb_filepath: str):
    """Extrai informações gerais do cabeçalho de um arquivo PDB (título, organismo, etc.)."""
    info = defaultdict(str) # Usando defaultdict para facilitar a concatenação de strings.
    try:
        with open(pdb_filepath, 'r') as f:
            for line in f:
                if line.startswith("HEADER"):
                    info['classification'] = line[10:50].strip()
                    info['dep_date'] = line[50:59].strip()
                elif line.startswith("TITLE"):
                    # O título pode ocupar várias linhas, então eu concateno.
                    info['title'] += line[10:80].strip() + " "
                elif line.startswith("SOURCE"):
                    # Procuro especificamente pela linha com o nome científico do organismo.
                    if "ORGANISM_SCIENTIFIC" in line:
                        info['organism'] += line.split(":")[-1].strip().replace(';','') + " "
    except Exception as e:
        print(f"Aviso: Não foi possível extrair informações do cabeçalho do PDB: {e}", file=sys.stderr)
    
    # Limpo espaços extras que podem ter sido criados pela concatenação.
    for key in info:
        info[key] = ' '.join(info[key].split())
    return info
    
def write_csv(filepath, header, data_rows):
    """Função auxiliar para escrever dados em um arquivo CSV."""
    try:
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(header) # Escreve o cabeçalho.
            writer.writerows(data_rows) # Escreve todas as linhas de dados.
        print(f"Resultados salvos com sucesso em '{filepath}'", file=sys.stderr)
    except IOError as e:
        print(f"Erro ao salvar o arquivo: {e}", file=sys.stderr)

def generate_sphere_points(n_points: int):
    """Gera pontos distribuídos uniformemente na superfície de uma esfera (algoritmo Fibonacci sphere)."""
    # Essa função é crucial para o cálculo do SASA.
    points = []
    if n_points <= 0: return points
    phi = (1 + math.sqrt(5)) / 2  # Proporção áurea.
    for i in range(n_points):
        y = 1 - (2 * i / (n_points - 1)) if n_points > 1 else 0
        radius = math.sqrt(1 - y * y)
        theta = 2 * math.pi * i / phi
        points.append((math.cos(theta) * radius, y, math.sin(theta) * radius))
    return points

# Funções de Download, Conversão e Análise (o coração da ferramenta)

def filter_pdb_content(pdb_filepath, chains=None, protein_only=False):
    """
    Filtra o conteúdo de um arquivo PDB baseado em critérios especificados.
    
    Args:
        pdb_filepath: Caminho para o arquivo PDB
        chains: Lista de IDs de cadeias a manter (ex: ['A', 'B']). Se None, mantém todas.
        protein_only: Se True, remove água (HOH), ligantes e heteroátomos (mantém apenas ATOM)
    
    Returns:
        Número de átomos mantidos após filtragem
    """
    filtered_lines = []
    atoms_kept = 0
    atoms_removed = 0
    
    with open(pdb_filepath, 'r') as f:
        for line in f:
            keep_line = True
            
            # Processa linhas ATOM e HETATM
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Se protein_only, mantém apenas linhas ATOM
                if protein_only and line.startswith("HETATM"):
                    keep_line = False
                    atoms_removed += 1
                
                # Se protein_only, remove água (HOH)
                if protein_only and line[17:20].strip() == "HOH":
                    keep_line = False
                    atoms_removed += 1
                
                # Filtra por chains se especificado
                if keep_line and chains is not None:
                    chain_id = line[21]
                    if chain_id not in chains:
                        keep_line = False
                        atoms_removed += 1
                
                if keep_line:
                    atoms_kept += 1
            
            # Mantém todas as outras linhas (HEADER, TITLE, etc.) exceto HETATM se protein_only
            elif not (protein_only and line.startswith("HETATM")):
                pass  # Mantém a linha
            else:
                keep_line = False
            
            if keep_line:
                filtered_lines.append(line)
    
    # Reescreve o arquivo com conteúdo filtrado
    with open(pdb_filepath, 'w') as f:
        f.writelines(filtered_lines)
    
    return atoms_kept, atoms_removed

def handle_fetch_pdb(args):
    """Baixa um arquivo PDB do banco de dados RCSB PDB."""
    pdb_id = args.pdb_id.upper()
    # Validação simples do ID do PDB.
    if len(pdb_id) != 4:
        print(f"Erro: O ID do PDB '{pdb_id}' é inválido. Deve conter 4 caracteres.", file=sys.stderr)
        return
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_file = args.output if args.output else f"{pdb_id}.pdb"
    try:
        print(f"Baixando {pdb_id} de {url}...", file=sys.stderr)
        # Faço o download e salvo o arquivo.
        with urllib.request.urlopen(url) as response, open(output_file, 'wb') as out_file:
            data = response.read()
            if not data: raise ValueError("Arquivo PDB recebido está vazio.") # Checa se o download não falhou silenciosamente.
            out_file.write(data)
        print(f"Arquivo PDB salvo com sucesso em '{output_file}'", file=sys.stderr)
        
        # Aplica filtros se especificados
        chains_list = None
        if args.chains:
            chains_list = [c.strip().upper() for c in args.chains.split(',')]
        
        if chains_list or args.protein_only:
            print("\n--- Aplicando Filtros ---", file=sys.stderr)
            if chains_list:
                print(f"  Cadeias selecionadas: {', '.join(chains_list)}", file=sys.stderr)
            if args.protein_only:
                print(f"  Modo: Apenas proteína (removendo água e ligantes)", file=sys.stderr)
            
            atoms_kept, atoms_removed = filter_pdb_content(output_file, chains_list, args.protein_only)
            print(f"  ✓ Filtragem concluída: {atoms_kept} átomos mantidos, {atoms_removed} removidos", file=sys.stderr)
        
        # Após o download, exibo algumas informações úteis do cabeçalho.
        info = extract_pdb_header_info(output_file)
        if info:
            print("\n--- Informações da Estrutura ---", file=sys.stderr)
            if info['title']: print(f"  Título    : {info['title']}", file=sys.stderr)
            if info['dep_date']: print(f"  Data      : {info['dep_date']}", file=sys.stderr)
            if info['classification']: print(f"  Classe    : {info['classification']}", file=sys.stderr)
            if info['organism']: print(f"  Organismo : {info['organism']}", file=sys.stderr)

    except urllib.error.HTTPError as e:
        print(f"Erro ao baixar o arquivo: Não foi possível encontrar o PDB ID '{pdb_id}'. (HTTP {e.code})", file=sys.stderr)
        if os.path.exists(output_file): os.remove(output_file) # Limpa o arquivo vazio se o download falhar.
    except Exception as e:
        print(f"Ocorreu um erro inesperado: {e}", file=sys.stderr)
        if os.path.exists(output_file): os.remove(output_file)

def handle_pdb_to_fasta(args):
    """Converte um arquivo PDB para o formato FASTA."""
    sequence = get_sequence_from_pdb(args.pdb_file)
    if sequence:
        # Monto o cabeçalho FASTA com o nome do arquivo de origem.
        header = f">sequence_from_{os.path.basename(args.pdb_file)}"
        fasta_output = f"{header}\n{sequence}\n"
        # Se o usuário especificou um arquivo de saída, salvo nele. Senão, imprimo na tela.
        if args.output:
            with open(args.output, 'w') as f: f.write(fasta_output)
            print(f"Sequência FASTA salva em '{args.output}'", file=sys.stderr)
        else:
            print(fasta_output)

def handle_csv_to_fasta(args):
    """Converte um arquivo CSV com IDs e sequências para o formato FASTA."""
    try:
        with open(args.csv_file, mode='r', newline='', encoding='utf-8') as infile:
            reader = csv.reader(infile, delimiter=args.delimiter)
            id_col_idx, seq_col_idx = -1, -1
            # Se o arquivo tem cabeçalho, procuro as colunas pelo nome ou índice.
            if args.header:
                header_row = next(reader)
                try: # Tenta converter para inteiro (índice)
                    id_col_idx, seq_col_idx = int(args.id_col), int(args.seq_col)
                except ValueError: # Se não for inteiro, é um nome de coluna.
                    if args.id_col in header_row: id_col_idx = header_row.index(args.id_col)
                    if args.seq_col in header_row: seq_col_idx = header_row.index(args.seq_col)
                if id_col_idx == -1 or seq_col_idx == -1:
                    print(f"Erro: Coluna de ID ('{args.id_col}') ou Sequência ('{args.seq_col}') não encontrada.", file=sys.stderr)
                    return
            else: # Se não tem cabeçalho, as colunas devem ser fornecidas como índices numéricos.
                try:
                    id_col_idx, seq_col_idx = int(args.id_col), int(args.seq_col)
                except ValueError:
                    print("Erro: Sem cabeçalho, as colunas devem ser índices numéricos.", file=sys.stderr)
                    return
            
            # Gero as strings FASTA para cada linha do CSV.
            fasta_records = [f">{row[id_col_idx].strip()}\n{row[seq_col_idx].strip().replace(' ', '')}" for row in reader if row and len(row) > max(id_col_idx, seq_col_idx)]
            output_content = "\n".join(fasta_records)
            
            # Salvo em arquivo ou imprimo na tela.
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
    """Calcula várias propriedades físico-químicas de uma sequência de proteína."""
    sequence = args.sequence.upper()
    if not sequence:
        print("Erro: Sequência de entrada está vazia.", file=sys.stderr)
        return
    
    length = len(sequence)
    aa_composition = {aa: sequence.count(aa) for aa in MOLECULAR_WEIGHT.keys()}
    
    # Peso Molecular: somo os pesos dos aminoácidos e subtraio o peso das moléculas de água perdidas na ligação peptídica.
    mw = sum(MOLECULAR_WEIGHT.get(aa, 0) * count for aa, count in aa_composition.items()) - (length - 1) * 18.015
    
    # GRAVY: média dos valores de hidropaticidade de Kyte-Doolittle.
    gravy = sum(KYTE_DOOLITTLE.get(aa, 0) for aa in sequence) / length if length > 0 else 0
    
    # Ponto Isoelétrico (pI): procuro o pH onde a carga líquida da proteína é mais próxima de zero.
    pi = 0.0
    min_charge_diff = float('inf')
    for ph_int in range(1401): # Itera de pH 0 a 14.00 com passo de 0.01.
        ph = ph_int * 0.01
        # Carga dos terminais N e C.
        net_charge = (10**PKA_VALUES['N-term']) / (10**PKA_VALUES['N-term'] + 10**ph) - (10**ph) / (10**PKA_VALUES['C-term'] + 10**ph)
        # Carga das cadeias laterais ionizáveis.
        for aa, count in aa_composition.items():
            if aa in ['R', 'H', 'K']: net_charge += count * (10**PKA_VALUES[aa]) / (10**PKA_VALUES[aa] + 10**ph)
            elif aa in ['D', 'E', 'C', 'Y']: net_charge -= count * (10**ph) / (10**PKA_VALUES[aa] + 10**ph)
        
        if abs(net_charge) < min_charge_diff:
            min_charge_diff, pi = abs(net_charge), ph
            
    # Índice Alifático: volume relativo ocupado por cadeias laterais alifáticas.
    a, b = 2.9, 3.9
    aliphatic_index = (aa_composition.get('A', 0) + a * aa_composition.get('V', 0) + b * (aa_composition.get('I', 0) + aa_composition.get('L', 0))) / length * 100 if length > 0 else 0
    
    # Índice de Instabilidade: uma estimativa da estabilidade da proteína in vivo.
    instability_index = (10 / length) * sum(DIWV[sequence[i]][sequence[i+1]] for i in range(length - 1)) if length > 1 else 0
    stability = "Estável" if instability_index < 40 else "Instável"
    
    # Meia-vida: regra do N-terminal. Uma estimativa grosseira.
    n_term = sequence[0]
    half_life_mammal = ">10 horas" if n_term in ['A','C','G','M','P','S','T','V'] else "2-30 min" if n_term in ['I','L','F','W','Y','D','E','N','Q'] else "Desconhecido"
    
    # Contagem de tipos de resíduos.
    def count_group(group_aas): return sum(aa_composition.get(aa, 0) for aa in group_aas)
    acidic, basic = count_group({'D', 'E'}), count_group({'R', 'K', 'H'})
    polar, non_polar = count_group({'N', 'Q', 'S', 'T', 'Y', 'C'}), count_group({'A', 'V', 'L', 'I', 'P', 'F', 'W', 'M', 'G'})
    
    # Formatação da saída.
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

    # Gera visualizações se solicitado
    if hasattr(args, 'plot_treemap') and args.plot_treemap:
        if HAS_VIZ:
            biohub_viz.plot_aa_composition_treemap(aa_composition, length, args.plot_treemap)
        else:
            print("Aviso: Módulo de visualização não disponível. Instale matplotlib, numpy e squarify.", file=sys.stderr)

    if hasattr(args, 'plot_composition') and args.plot_composition:
        if HAS_VIZ:
            biohub_viz.plot_aa_composition_bar(aa_composition, length, args.plot_composition)
        else:
            print("Aviso: Módulo de visualização não disponível. Instale matplotlib.", file=sys.stderr)

    if hasattr(args, 'plot_hydro') and args.plot_hydro:
        if HAS_VIZ:
            window = args.window if hasattr(args, 'window') else 9
            biohub_viz.plot_hydropathy_profile(sequence, KYTE_DOOLITTLE, args.plot_hydro, window)
        else:
            print("Aviso: Módulo de visualização não disponível. Instale matplotlib.", file=sys.stderr)

def calculate_intramolecular_contacts(args):
    """Calcula contatos entre resíduos com base na distância mínima entre quaisquer átomos."""
    # Dicionário: resíduo -> lista de coordenadas de todos os átomos
    residue_atoms = {}
    if not os.path.exists(args.pdb_file):
        print(f"Erro: Arquivo não encontrado em '{args.pdb_file}'", file=sys.stderr)
        return
    with open(args.pdb_file, 'r') as f:
        first_chain_id = None
        for line in f:
            # Pego todos os átomos (não apenas CA) da primeira cadeia
            if line.startswith("ATOM"):
                if first_chain_id is None: first_chain_id = line[21]
                if line[21] == first_chain_id:
                    res_num = int(line[22:26])
                    coords = tuple(float(line[i:i+8]) for i in [30, 38, 46])
                    if res_num not in residue_atoms:
                        residue_atoms[res_num] = []
                    residue_atoms[res_num].append(coords)

    residues = sorted(residue_atoms.keys())

    # Calcula a distância mínima entre qualquer par de átomos de dois resíduos
    contacts = []
    for i, r1 in enumerate(residues):
        for r2 in residues[i+1:]:
            # Ignora vizinhos diretos na sequência
            if abs(r1 - r2) <= 1:
                continue

            # Calcula distância mínima entre qualquer átomo de r1 e qualquer átomo de r2
            min_dist = float('inf')
            for atom1 in residue_atoms[r1]:
                for atom2 in residue_atoms[r2]:
                    dist = math.sqrt(sum((c1-c2)**2 for c1,c2 in zip(atom1, atom2)))
                    if dist < min_dist:
                        min_dist = dist

            # Se a distância mínima está abaixo do threshold, é um contato
            if min_dist <= args.threshold:
                contacts.append((r1, r2, round(min_dist, 3)))
    
    if args.output:
        write_csv(args.output, ["Residuo_1", "Residuo_2", "Distancia_A"], contacts)
    else:
        print(f"--- Contatos Intramoleculares (Limiar = {args.threshold:.1f} Å) ---")
        if not contacts: print("Nenhum contato encontrado.")
        else:
            for c in contacts: print(f"Res {c[0]} - Res {c[1]}: {c[2]:.3f} Å")

    # Gera visualização se solicitado
    if hasattr(args, 'plot') and args.plot and contacts:
        if HAS_VIZ:
            max_res = max(max(c[0], c[1]) for c in contacts) if contacts else 1
            biohub_viz.plot_contact_map(contacts, max_res, args.plot, args.threshold)
        else:
            print("Aviso: Módulo de visualização não disponível. Instale matplotlib e numpy.", file=sys.stderr)

def write_pdb_with_bfactor(input_pdb_path, output_pdb_path, atom_values, property_name="Property"):
    """
    Reescreve um arquivo PDB substituindo os valores do B-factor por valores calculados.
    
    Args:
        input_pdb_path: Caminho do PDB original
        output_pdb_path: Caminho para salvar o PDB anotado
        atom_values: Lista de dicionários com 'atom_num' e 'value'
        property_name: Nome da propriedade sendo escrita (para mensagens)
    """
    # Cria um dicionário para lookup rápido por número de átomo
    value_by_atom = {atom['atom_num']: atom['value'] for atom in atom_values}
    
    # Calcula min/max para estatísticas
    values_list = [atom['value'] for atom in atom_values]
    min_val = min(values_list)
    max_val = max(values_list)
    
    atoms_updated = 0
    
    try:
        with open(input_pdb_path, 'r') as infile, open(output_pdb_path, 'w') as outfile:
            for line in infile:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    # Extrai o número do átomo
                    atom_num = int(line[6:11].strip())
                    
                    if atom_num in value_by_atom:
                        # Pega o valor calculado
                        new_bfactor = value_by_atom[atom_num]
                        
                        # Reconstrói a linha com novo B-factor (colunas 61-66)
                        # Formato: %6.2f (6 caracteres, 2 decimais)
                        new_line = line[:60] + f"{new_bfactor:6.2f}" + line[66:]
                        outfile.write(new_line)
                        atoms_updated += 1
                    else:
                        # Átomo não encontrado nos valores calculados, mantém original
                        outfile.write(line)
                else:
                    # Linhas que não são ATOM/HETATM, mantém inalteradas
                    outfile.write(line)
        
        print(f"PDB anotado salvo em '{output_pdb_path}'", file=sys.stderr)
        print(f"  {atoms_updated} átomos tiveram o B-factor atualizado com {property_name}", file=sys.stderr)
        print(f"  Range de valores: {min_val:.2f} - {max_val:.2f}", file=sys.stderr)
        
    except Exception as e:
        print(f"Erro ao escrever PDB anotado: {e}", file=sys.stderr)

def generate_pymol_session(pdb_path, output_pse, property_type="hydrophobicity", min_val=-4.5, max_val=4.5):
    """
    Gera um arquivo de sessão PyMOL (.pse) com visualização configurada.
    Se PyMOL não estiver disponível, gera um script .pml que pode ser executado manualmente.
    
    Args:
        pdb_path: Caminho do PDB anotado (com B-factor modificado)
        output_pse: Caminho para salvar o arquivo .pse
        property_type: Tipo de propriedade ('hydrophobicity' ou 'sasa')
        min_val: Valor mínimo para o gradiente de cores
        max_val: Valor máximo para o gradiente de cores
    """
    try:
        # Gera também um script .pml independente do resultado
        pml_script = output_pse.replace('.pse', '.pml')
        
        # Cria um script PyMOL
        script_content = []
        
        # Carrega o PDB
        pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
        script_content.append(f"# Script PyMOL gerado pelo BioHub")
        script_content.append(f"# Propriedade: {property_type}")
        script_content.append(f"")
        script_content.append(f"# Carrega a estrutura")
        script_content.append(f"load {os.path.abspath(pdb_path)}, {pdb_name}")
        script_content.append(f"")
        script_content.append(f"# Remove todas as representações padrão")
        script_content.append(f"hide everything, {pdb_name}")
        script_content.append(f"")
        
        # Configurações de visualização baseadas no tipo de propriedade
        if property_type == "hydrophobicity":
            # Esquema de cores: azul (hidrofílico) -> branco -> vermelho (hidrofóbico)
            script_content.append(f"# === HIDROFOBICIDADE ===")
            script_content.append(f"# Azul = Hidrofílico ({min_val}), Vermelho = Hidrofóbico ({max_val})")
            script_content.append(f"")
            script_content.append(f"# Representação Cartoon (fita)")
            script_content.append(f"show cartoon, {pdb_name}")
            script_content.append(f"cartoon automatic, {pdb_name}")
            script_content.append(f"set cartoon_fancy_helices, 1")
            script_content.append(f"spectrum b, blue_white_red, {pdb_name}, minimum={min_val}, maximum={max_val}")
            script_content.append(f"")
            script_content.append(f"# Representação Sticks (bastões)")
            script_content.append(f"show sticks, {pdb_name}")
            script_content.append(f"set stick_radius, 0.2, {pdb_name}")
            script_content.append(f"set stick_color, gray, {pdb_name}")
            script_content.append(f"util.cbag {pdb_name}")  # Cores por átomo (C=cinza, N=azul, O=vermelho)
            script_content.append(f"")
            script_content.append(f"# Representação Surface (superfície)")
            script_content.append(f"show surface, {pdb_name}")
            script_content.append(f"set surface_quality, 1")
            script_content.append(f"set transparency, 0.5, {pdb_name}")
            script_content.append(f"# Aplica gradiente de hidrofobicidade na superfície")
            script_content.append(f"set surface_color, white, {pdb_name}")
            script_content.append(f"spectrum b, blue_white_red, {pdb_name}, minimum={min_val}, maximum={max_val}")
            
        elif property_type == "sasa":
            # Esquema de cores INVERTIDO: vermelho (enterrado) -> branco -> azul (exposto)
            # Invertemos porque alto SASA = exposto ao solvente (água) = deve ser azul
            script_content.append(f"# === SASA (Acessibilidade ao Solvente) ===")
            script_content.append(f"# Vermelho = Enterrado ({min_val:.1f} Ų), Azul = Exposto ({max_val:.1f} Ų)")
            script_content.append(f"")
            script_content.append(f"# Representação Cartoon (fita)")
            script_content.append(f"show cartoon, {pdb_name}")
            script_content.append(f"cartoon automatic, {pdb_name}")
            script_content.append(f"set cartoon_fancy_helices, 1")
            script_content.append(f"spectrum b, red_white_blue, {pdb_name}, minimum={min_val}, maximum={max_val}")
            script_content.append(f"")
            script_content.append(f"# Representação Sticks (bastões)")
            script_content.append(f"show sticks, {pdb_name}")
            script_content.append(f"set stick_radius, 0.2, {pdb_name}")
            script_content.append(f"set stick_color, gray, {pdb_name}")
            script_content.append(f"util.cbag {pdb_name}")  # Cores por átomo
            script_content.append(f"")
            script_content.append(f"# Representação Surface (superfície)")
            script_content.append(f"show surface, {pdb_name}")
            script_content.append(f"set surface_quality, 1")
            script_content.append(f"set transparency, 0.5, {pdb_name}")
            script_content.append(f"# Aplica gradiente de SASA na superfície")
            script_content.append(f"set surface_color, white, {pdb_name}")
            script_content.append(f"spectrum b, red_white_blue, {pdb_name}, minimum={min_val}, maximum={max_val}")
        
        script_content.append(f"")
        script_content.append(f"# === Configurações Gerais de Qualidade ===")
        script_content.append(f"bg_color white")
        script_content.append(f"set ray_shadow, 0")
        script_content.append(f"set antialias, 2")
        script_content.append(f"set orthoscopic, 0")
        script_content.append(f"set valence, 0")
        script_content.append(f"")
        script_content.append(f"# Centra e ajusta visualização")
        script_content.append(f"center {pdb_name}")
        script_content.append(f"zoom {pdb_name}")
        script_content.append(f"orient {pdb_name}")
        script_content.append(f"")
        script_content.append(f"# Comandos úteis:")
        script_content.append(f"# hide surface - ocultar superfície")
        script_content.append(f"# hide sticks - ocultar bastões")
        script_content.append(f"# hide cartoon - ocultar cartoon")
        script_content.append(f"# set transparency, 0.7 - aumentar transparência")
        script_content.append(f"")
        script_content.append(f"# Para salvar a sessão, use:")
        script_content.append(f"# save {os.path.abspath(output_pse)}")
        
        # Salva o script .pml
        with open(pml_script, 'w') as f:
            f.write('\n'.join(script_content))
        
        print(f"Script PyMOL salvo em '{pml_script}'", file=sys.stderr)
        print(f"  Execute com: pymol {pml_script}", file=sys.stderr)
        
        # Tenta gerar o .pse automaticamente se PyMOL estiver disponível
        if shutil.which('pymol'):
            # Adiciona comandos para salvar sessão e sair
            save_script_content = script_content + [
                f"",
                f"save {os.path.abspath(output_pse)}",
                f"quit"
            ]
            
            # Escreve script temporário
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pml', delete=False) as f:
                script_path = f.name
                f.write('\n'.join(save_script_content))
            
            try:
                result = subprocess.run(
                    ['pymol', '-c', '-q', script_path],
                    capture_output=True,
                    text=True,
                    timeout=30
                )
                
                if result.returncode == 0:
                    print(f"Sessão PyMOL salva em '{output_pse}'", file=sys.stderr)
                    print(f"  Abra com: pymol {output_pse}", file=sys.stderr)
                else:
                    print(f"Aviso: Não foi possível gerar o arquivo .pse automaticamente.", file=sys.stderr)
                    print(f"  Use o script .pml manualmente: pymol {pml_script}", file=sys.stderr)
            
            except subprocess.TimeoutExpired:
                print(f"Aviso: Timeout ao executar PyMOL.", file=sys.stderr)
            finally:
                # Remove script temporário
                if os.path.exists(script_path):
                    os.remove(script_path)
        else:
            print(f"  Nota: PyMOL não está no PATH. Instale-o para gerar o arquivo .pse automaticamente.", file=sys.stderr)
        
    except Exception as e:
        print(f"Erro ao gerar script PyMOL: {e}", file=sys.stderr)

def predict_solvent_hydrophoby(args):
    """Prevê a exposição ao solvente usando hidrofobicidade (Kyte-Doolittle) por átomo."""
    atoms = parse_pdb_atoms(args.pdb_file)
    if not atoms: return
    
    results = []
    
    # Para cada átomo, atribui a hidrofobicidade do seu resíduo
    for atom in atoms:
        res_name = atom["res_name"]
        # Converte nome de 3 letras para 1 letra
        aa_code = THREE_TO_ONE.get(res_name, None)
        
        if aa_code:
            # Pega o score de hidrofobicidade do aminoácido
            hydro_score = KYTE_DOOLITTLE.get(aa_code, 0.0)
        else:
            # Se não for aminoácido padrão (ex: ligante), define como 0
            hydro_score = 0.0
        
        results.append({
            "chain_id": atom["chain_id"],
            "res_num": atom["res_num"],
            "res_name": atom["res_name"],
            "atom_num": atom["atom_num"],
            "atom_name": atom["atom_name"],
            "hydrophobicity": hydro_score
        })
    
    print(f"Total de átomos analisados: {len(results)}", file=sys.stderr)
    
    # Gera dados para saída
    results_data = [
        [atom["chain_id"], atom["res_num"], atom["res_name"], atom["atom_num"], atom["atom_name"], f"{atom['hydrophobicity']:.3f}"]
        for atom in results
    ]
    
    if args.output:
        write_csv(args.output, ["Chain", "ResNum", "ResName", "AtomNum", "AtomName", "Hydrophobicity"], results_data)
    else:
        print(f"--- Hidrofobicidade por Átomo (Escala Kyte-Doolittle) ---")
        print("Chain | ResNum | ResName | AtomNum | AtomName | Hydrophobicity")
        for row in results_data[:20]:  # Mostra apenas os primeiros 20
            print(f"{row[0]:<5} | {row[1]:<6} | {row[2]:<7} | {row[3]:<7} | {row[4]:<8} | {row[5]}")
        if len(results_data) > 20:
            print(f"... e mais {len(results_data) - 20} átomos. Use -o para salvar todos os dados.")
    
    # Gera visualização se solicitado
    if hasattr(args, 'plot_hydrophoby') and args.plot_hydrophoby:
        if HAS_VIZ:
            biohub_viz.plot_hydrophoby_profile(results, args.plot_hydrophoby)
        else:
            print("Aviso: Módulo de visualização não disponível. Instale matplotlib e numpy.", file=sys.stderr)
    
    # Gera PDB anotado se solicitado
    if args.write_pdb:
        # Prepara dados para escrita no B-factor
        atom_values = [{'atom_num': atom['atom_num'], 'value': atom['hydrophobicity']} for atom in results]
        write_pdb_with_bfactor(args.pdb_file, args.write_pdb, atom_values, "Hydrophobicity")
        
        # Gera sessão PyMOL se solicitado
        if args.pymol:
            generate_pymol_session(args.write_pdb, args.pymol, property_type="hydrophobicity", min_val=-4.5, max_val=4.5)
    elif args.pymol:
        # Se --pymol foi especificado mas --write-pdb não, avisa o usuário
        print("Aviso: --pymol requer --write-pdb. Gerando PDB temporário...", file=sys.stderr)
        temp_pdb = "temp_hydro.pdb"
        atom_values = [{'atom_num': atom['atom_num'], 'value': atom['hydrophobicity']} for atom in results]
        write_pdb_with_bfactor(args.pdb_file, temp_pdb, atom_values, "Hydrophobicity")
        generate_pymol_session(temp_pdb, args.pymol, property_type="hydrophobicity", min_val=-4.5, max_val=4.5)

def calculate_sasa(args):
    """Calcula a Área de Superfície Acessível ao Solvente (SASA) usando o método de Shrake-Rupley."""
    atoms = parse_pdb_atoms(args.pdb_file)
    if not atoms: return
    # Gero os pontos na esfera que serão usados para testar a acessibilidade de cada átomo.
    sphere_points = generate_sphere_points(args.num_points)
    total_sasa = 0.0
    sasa_per_atom = []  # Lista para armazenar SASA de cada átomo
    
    # Para cada átomo...
    for i, atom_i in enumerate(atoms):
        radius_i = VDW_RADII.get(atom_i["element"], VDW_RADII['DEFAULT'])
        extended_radius = radius_i + args.probe_radius # Raio do átomo + raio da sonda (água).
        accessible_points = 0
        # ... testo cada ponto na sua superfície estendida.
        for sp in sphere_points:
            point_is_accessible = True
            point = (atom_i["x"] + extended_radius * sp[0], atom_i["y"] + extended_radius * sp[1], atom_i["z"] + extended_radius * sp[2])
            # Verifico se este ponto está dentro da esfera estendida de qualquer outro átomo.
            for j, atom_j in enumerate(atoms):
                if i == j: continue
                radius_j = VDW_RADII.get(atom_j["element"], VDW_RADII['DEFAULT'])
                # Se a distância do ponto ao centro do átomo j for menor que o raio estendido de j, o ponto está ocluído.
                if sum((p-c)**2 for p,c in zip(point, (atom_j["x"], atom_j["y"], atom_j["z"]))) < (radius_j + args.probe_radius)**2:
                    point_is_accessible = False; break
            if point_is_accessible: accessible_points += 1
        
        # O SASA do átomo é a proporção de pontos acessíveis multiplicada pela área da esfera estendida.
        atom_sasa = (accessible_points / args.num_points) * 4.0 * math.pi * extended_radius**2 if args.num_points > 0 else 0
        total_sasa += atom_sasa
        
        # Armazena dados do átomo com SASA
        sasa_per_atom.append({
            "chain_id": atom_i["chain_id"],
            "res_num": atom_i["res_num"],
            "res_name": atom_i["res_name"],
            "atom_num": atom_i["atom_num"],
            "atom_name": atom_i["atom_name"],
            "sasa": atom_sasa
        })
        
    print(f"SASA Total da Molécula: {total_sasa:.2f} Å²", file=sys.stderr)
    print(f"Total de átomos analisados: {len(sasa_per_atom)}", file=sys.stderr)
    
    # Gera dados para saída
    results_data = [
        [atom["chain_id"], atom["res_num"], atom["res_name"], atom["atom_num"], atom["atom_name"], f"{atom['sasa']:.2f}"]
        for atom in sasa_per_atom
    ]
    
    if args.output:
        write_csv(args.output, ["Chain", "ResNum", "ResName", "AtomNum", "AtomName", "SASA_A2"], results_data)
    else:
        print("--- SASA por Átomo ---")
        print("Chain | ResNum | ResName | AtomNum | AtomName | SASA (Å²)")
        for row in results_data[:20]:  # Mostra apenas os primeiros 20 para não poluir o terminal
            print(f"{row[0]:<5} | {row[1]:<6} | {row[2]:<7} | {row[3]:<7} | {row[4]:<8} | {row[5]}")
        if len(results_data) > 20:
            print(f"... e mais {len(results_data) - 20} átomos. Use -o para salvar todos os dados.")

    # Gera visualização se solicitado
    if hasattr(args, 'plot_profile') and args.plot_profile:
        if HAS_VIZ:
            biohub_viz.plot_sasa_profile(sasa_per_atom, args.plot_profile)
        else:
            print("Aviso: Módulo de visualização não disponível. Instale matplotlib e numpy.", file=sys.stderr)

    # Gera PDB anotado se solicitado
    if args.write_pdb:
        # Calcula SASA MÉDIO POR RESÍDUO para visualização mais biologicamente relevante
        from collections import defaultdict
        residue_sasa = defaultdict(list)
        
        # Agrupa SASA por resíduo (chain + res_num)
        for atom in sasa_per_atom:
            residue_key = (atom['chain_id'], atom['res_num'])
            residue_sasa[residue_key].append(atom['sasa'])
        
        # Calcula média por resíduo
        residue_avg_sasa = {}
        for residue_key, sasa_list in residue_sasa.items():
            residue_avg_sasa[residue_key] = sum(sasa_list) / len(sasa_list)
        
        # Atribui a média do resíduo a todos os átomos daquele resíduo
        atom_values = []
        for atom in sasa_per_atom:
            residue_key = (atom['chain_id'], atom['res_num'])
            avg_sasa = residue_avg_sasa[residue_key]
            atom_values.append({'atom_num': atom['atom_num'], 'value': avg_sasa})
        
        write_pdb_with_bfactor(args.pdb_file, args.write_pdb, atom_values, "SASA (média por resíduo)")
        
        # Gera sessão PyMOL se solicitado
        if args.pymol:
            # Para SASA, usa a média por resíduo para definir o range
            avg_sasa_values = sorted(list(residue_avg_sasa.values()))
            
            # Remove valores muito baixos para análise (resíduos quase completamente enterrados)
            non_zero_sasa = [v for v in avg_sasa_values if v > 0.5]
            
            if len(non_zero_sasa) > 0:
                # Usa percentil 70 dos valores não-zero para melhor sensibilidade
                percentil_70_idx = int(len(non_zero_sasa) * 0.70)
                max_sasa = non_zero_sasa[percentil_70_idx]
                min_sasa = 0.0
            else:
                min_sasa = 0.0
                max_sasa = max(avg_sasa_values) if avg_sasa_values else 1.0
            
            print(f"Range de visualização SASA (média por resíduo): 0.00 - {max_sasa:.2f} Ų (percentil 70)", file=sys.stderr)
            print(f"  Resíduos totais: {len(residue_avg_sasa)}", file=sys.stderr)
            print(f"  Resíduos enterrados (SASA<0.5): {len(avg_sasa_values) - len(non_zero_sasa)}", file=sys.stderr)
            print(f"  Resíduos expostos (SASA≥0.5): {len(non_zero_sasa)}", file=sys.stderr)
            generate_pymol_session(args.write_pdb, args.pymol, property_type="sasa", min_val=min_sasa, max_val=max_sasa)
    elif args.pymol:
        # Se --pymol foi especificado mas --write-pdb não, avisa o usuário
        print("Aviso: --pymol requer --write-pdb. Gerando PDB temporário...", file=sys.stderr)
        
        # Calcula média por resíduo também para o caso temporário
        from collections import defaultdict
        residue_sasa = defaultdict(list)
        for atom in sasa_per_atom:
            residue_key = (atom['chain_id'], atom['res_num'])
            residue_sasa[residue_key].append(atom['sasa'])
        
        residue_avg_sasa = {}
        for residue_key, sasa_list in residue_sasa.items():
            residue_avg_sasa[residue_key] = sum(sasa_list) / len(sasa_list)
        
        atom_values = []
        for atom in sasa_per_atom:
            residue_key = (atom['chain_id'], atom['res_num'])
            avg_sasa = residue_avg_sasa[residue_key]
            atom_values.append({'atom_num': atom['atom_num'], 'value': avg_sasa})
        
        temp_pdb = "temp_sasa.pdb"
        write_pdb_with_bfactor(args.pdb_file, temp_pdb, atom_values, "SASA (média por resíduo)")
        
        avg_sasa_values = sorted(list(residue_avg_sasa.values()))
        non_zero_sasa = [v for v in avg_sasa_values if v > 0.5]
        
        if len(non_zero_sasa) > 0:
            percentil_70_idx = int(len(non_zero_sasa) * 0.70)
            max_sasa = non_zero_sasa[percentil_70_idx]
            min_sasa = 0.0
        else:
            min_sasa = 0.0
            max_sasa = max(avg_sasa_values) if avg_sasa_values else 1.0
        
        print(f"Range de visualização SASA (média por resíduo): 0.00 - {max_sasa:.2f} Ų (percentil 70)", file=sys.stderr)
        print(f"  Resíduos totais: {len(residue_avg_sasa)}", file=sys.stderr)
        print(f"  Resíduos enterrados (SASA<0.5): {len(avg_sasa_values) - len(non_zero_sasa)}", file=sys.stderr)
        print(f"  Resíduos expostos (SASA≥0.5): {len(non_zero_sasa)}", file=sys.stderr)
        generate_pymol_session(temp_pdb, args.pymol, property_type="sasa", min_val=min_sasa, max_val=max_sasa)

def run_apbs_analysis(args): #BETA, TALVEZ SERÁ DESCONTINUADO
    """Executa PDB2PQR e APBS para calcular a energia de solvatação eletrostática."""
    # Verifico se os programas externos necessários estão instalados e no PATH do sistema.
    for exe in ["pdb2pqr", "apbs"]:
        if not shutil.which(exe):
            print(f"Erro: '{exe}' não encontrado no PATH.", file=sys.stderr); return
            
    work_dir = tempfile.mkdtemp() # Crio um diretório temporário para os arquivos.
    try:
        # 1. Converter PDB para PQR usando pdb2pqr.
        pqr_path = os.path.join(work_dir, f"{os.path.basename(args.pdb_file)}.pqr")
        print(f"1. Executando pdb2pqr...", file=sys.stderr)
        subprocess.run(["pdb2pqr", "--ff=amber", args.pdb_file, pqr_path], check=True, capture_output=True, text=True, timeout=60)
        
        # 2. Criar o arquivo de entrada para o APBS.
        apbs_in_path = os.path.join(work_dir, "apbs.in")
        with open(apbs_in_path, 'w') as f: f.write(f"read\n    mol pqr {os.path.basename(pqr_path)}\nend\nelec\n    mg-auto\nend\nquit\n")
        
        # 3. Executar o APBS.
        print("2. Executando APBS...", file=sys.stderr)
        result = subprocess.run(["apbs", "apbs.in"], check=True, capture_output=True, text=True, cwd=work_dir, timeout=300)
        
        # 4. Extrair a energia de solvatação do resultado.
        print("3. Analisando resultados...", file=sys.stderr)
        energy = "Não encontrada"
        for line in result.stdout.splitlines():
            if "Global net ELEC energy" in line:
                energy = f"{float(line.split()[-2]):.4f} kJ/mol"; break
                
        print(f"--- Energia de Solvatação Eletrostática ---")
        print(f"Energia: {energy}")
        
    except subprocess.CalledProcessError as e:
        # Se algum dos programas externos falhar, mostro o erro.
        print(f"Erro ao executar processo externo:\n{e.stderr}", file=sys.stderr)
    finally:
        # Limpo os arquivos temporários, a menos que o usuário peça para mantê-los.
        if not args.no_cleanup: shutil.rmtree(work_dir)
        else: print(f"Arquivos intermediários mantidos em: {work_dir}", file=sys.stderr)

# Configuração da Interface de Linha de Comando

def main():
    """Função principal que organiza e executa a ferramenta."""
    # Exibo o banner se o script for executado interativamente ou sem argumentos.
    if sys.stdout.isatty() or len(sys.argv) == 1:
        print_banner()

    # Crio o parser principal com argparse.
    parser = argparse.ArgumentParser(
        description="BioHub: Uma ferramenta CLI para análise de proteínas.",
        formatter_class=argparse.RawTextHelpFormatter, # Para manter a formatação do texto de ajuda.
        epilog="Use: biohub.py [COMANDO] -h para ajuda detalhada sobre um comando."
    )
    # Crio os "subcomandos" (fetchpdb, fasta, etc.).
    subparsers = parser.add_subparsers(dest="command", help="\nComando a ser executado.")

    # Comando fetchpdb 
    parser_fetch = subparsers.add_parser("fetchpdb", help="Baixa um arquivo PDB do RCSB.", formatter_class=argparse.RawTextHelpFormatter)
    parser_fetch.add_argument("pdb_id", metavar="PDB_ID", help="O ID de 4 caracteres do PDB a ser baixado (ex: 1A2B).")
    parser_fetch.add_argument("-o", "--output", metavar="ARQUIVO", help="Nome do arquivo de saída (padrão: [PDB_ID].pdb).")
    parser_fetch.add_argument("--chains", metavar="CHAINS", help="Cadeias a serem mantidas, separadas por vírgula (ex: A,B). Se omitido, mantém todas.")
    parser_fetch.add_argument("--protein-only", action="store_true", help="Mantém apenas átomos de proteína (remove água, ligantes e heteroátomos).")

    # Comando fasta 
    parser_fasta = subparsers.add_parser("fasta", help="Converte um arquivo PDB em uma sequência FASTA.", formatter_class=argparse.RawTextHelpFormatter)
    parser_fasta.add_argument("pdb_file", metavar="ARQUIVO_PDB", help="Caminho para o arquivo PDB de entrada.")
    parser_fasta.add_argument("-o", "--output", metavar="ARQUIVO", help="Salva a saída em um arquivo FASTA (padrão: stdout).")
    
    # Comando csv2fasta 
    parser_csv = subparsers.add_parser("csv2fasta", help="Converte um arquivo CSV em um formato FASTA.", formatter_class=argparse.RawTextHelpFormatter)
    parser_csv.add_argument("csv_file", metavar="ARQUIVO_CSV", help="Caminho para o arquivo CSV de entrada.")
    parser_csv.add_argument("-o", "--output", metavar="ARQUIVO", help="Salva a saída em um arquivo FASTA (padrão: stdout).")
    parser_csv.add_argument("--id-col", metavar="COLUNA", default="0", help="Coluna do identificador (nome ou índice baseado em 0). Padrão: 0.")
    parser_csv.add_argument("--seq-col", metavar="COLUNA", default="1", help="Coluna da sequência (nome ou índice baseado em 0). Padrão: 1.")
    parser_csv.add_argument("--header", action="store_true", help="Flag para indicar que a primeira linha do CSV é um cabeçalho.")
    parser_csv.add_argument("--delimiter", metavar="CHAR", default=",", help="Caractere usado como delimitador no CSV. Padrão: ','.")

    # Comando physchem
    parser_physchem = subparsers.add_parser("physchem", help="Calcula um conjunto expandido de propriedades físico-químicas.", formatter_class=argparse.RawTextHelpFormatter)
    parser_physchem.add_argument("sequence", metavar="SEQUENCIA", help="A sequência de aminoácidos a ser analisada.")
    parser_physchem.add_argument("-o", "--output", metavar="ARQUIVO_CSV", help="Salva os resultados em um arquivo CSV.")
    parser_physchem.add_argument("--plot-treemap", metavar="ARQUIVO_PNG", help="Gera treemap de composição de aminoácidos (requer matplotlib, numpy, squarify).")
    parser_physchem.add_argument("--plot-composition", metavar="ARQUIVO_PNG", help="Gera gráfico de barras de composição de aminoácidos (requer matplotlib).")
    parser_physchem.add_argument("--plot-hydro", metavar="ARQUIVO_PNG", help="Gera perfil de hidrofobicidade Kyte-Doolittle (requer matplotlib).")
    parser_physchem.add_argument("--window", metavar="INT", type=int, default=9, help="Tamanho da janela para perfil de hidrofobicidade (padrão: 9).")

    # Comando contacts
    parser_contacts = subparsers.add_parser("contacts", help="Calcula contatos intramoleculares com base na distância entre C-Alfas.", formatter_class=argparse.RawTextHelpFormatter)
    parser_contacts.add_argument("pdb_file", metavar="ARQUIVO_PDB", help="Caminho para o arquivo PDB de entrada.")
    parser_contacts.add_argument("-t", "--threshold", metavar="FLOAT", type=float, default=8.0, help="Distância máxima em Angstroms para considerar um contato. Padrão: 8.0.")
    parser_contacts.add_argument("-o", "--output", metavar="ARQUIVO_CSV", help="Salva os resultados em um arquivo CSV.")
    parser_contacts.add_argument("--plot", metavar="ARQUIVO_PNG", help="Gera mapa de contatos (contact map) (requer matplotlib e numpy).")

    # Comando hydrophoby 
    parser_hydrophoby = subparsers.add_parser("hydrophoby", help="Calcula hidrofobicidade por átomo (escala Kyte-Doolittle).", formatter_class=argparse.RawTextHelpFormatter)
    parser_hydrophoby.add_argument("pdb_file", metavar="ARQUIVO_PDB", help="Caminho para o arquivo PDB de entrada.")
    parser_hydrophoby.add_argument("-o", "--output", metavar="ARQUIVO_CSV", help="Salva os resultados em um arquivo CSV.")
    parser_hydrophoby.add_argument("--write-pdb", metavar="ARQUIVO_PDB", help="Gera um arquivo PDB com a hidrofobicidade escrita no B-factor.")
    parser_hydrophoby.add_argument("--pymol", metavar="ARQUIVO_PSE", help="Gera um arquivo de sessão PyMOL (.pse) com visualização de hidrofobicidade.")
    parser_hydrophoby.add_argument("--plot-hydrophoby", metavar="ARQUIVO_PNG", help="Gera perfil de hidrofobicidade por resíduo (requer matplotlib e numpy).")
    
    # Comando sasa
    parser_sasa = subparsers.add_parser("sasa", help="Calcula a Área de Superfície Acessível ao Solvente (SASA).", formatter_class=argparse.RawTextHelpFormatter)
    parser_sasa.add_argument("pdb_file", metavar="ARQUIVO_PDB", help="Caminho para o arquivo PDB de entrada.")
    parser_sasa.add_argument("--probe-radius", metavar="FLOAT", type=float, default=1.4, help="Raio da sonda do solvente em Angstroms (padrão: 1.4 para água).")
    parser_sasa.add_argument("--num-points", metavar="INT", type=int, default=960, help="Número de pontos na superfície de cada átomo para o cálculo. Padrão: 960.")
    parser_sasa.add_argument("-o", "--output", metavar="ARQUIVO_CSV", help="Salva os resultados por átomo em um arquivo CSV.")
    parser_sasa.add_argument("--write-pdb", metavar="ARQUIVO_PDB", help="Gera um arquivo PDB com o SASA escrito no B-factor.")
    parser_sasa.add_argument("--pymol", metavar="ARQUIVO_PSE", help="Gera um arquivo de sessão PyMOL (.pse) com visualização de SASA.")
    parser_sasa.add_argument("--plot-profile", metavar="ARQUIVO_PNG", help="Gera perfil de SASA por resíduo (requer matplotlib e numpy).")

    # Comando apbs 
    parser_apbs = subparsers.add_parser("apbs", help="Calcula a energia de solvatação eletrostática (requer PDB2PQR e APBS).", formatter_class=argparse.RawTextHelpFormatter)
    parser_apbs.add_argument("pdb_file", metavar="ARQUIVO_PDB", help="Caminho para o arquivo PDB de entrada.")
    parser_apbs.add_argument("--no-cleanup", action="store_true", help="Previne a remoção dos arquivos temporários (PQR, apbs.in, etc.).")

    # Se nenhum comando for fornecido, exibo a ajuda e saio.
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    # Analiso os argumentos fornecidos pelo usuário.
    args = parser.parse_args()
    
    # Crio um dicionário para mapear o nome do comando à função que deve ser executada.
    command_functions = {
        "fetchpdb": handle_fetch_pdb, "fasta": handle_pdb_to_fasta, "csv2fasta": handle_csv_to_fasta,
        "physchem": calculate_physicochemical_properties, "contacts": calculate_intramolecular_contacts,
        "hydrophoby": predict_solvent_hydrophoby, "sasa": calculate_sasa, "apbs": run_apbs_analysis
    }
    
    # Chamo a função correspondente ao comando que o usuário escolheu.
    if args.command in command_functions:
        command_functions[args.command](args)
    else:
        # Se o comando não for válido, exibo a ajuda.
        parser.print_help(sys.stderr)

    # Se estiver em um terminal interativo, mostro a mensagem de saída.
    if sys.stdout.isatty():
        print_exit_message()

# Este é o ponto de entrada padrão para um script Python.
if __name__ == "__main__":
    main()
