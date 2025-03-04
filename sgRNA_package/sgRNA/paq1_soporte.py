import subprocess
from Bio import SeqIO
import pandas as pd
import re
from paq1_percent import predecir_eficiencia_nn,complemento_inverso, predecir_eficiencia, predecir_eficiencia_xgb

def load_file(file_path, file_type="fasta"):
    """Carga secuencias desde un archivo, compatible con formatos FASTA."""
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq).upper())
    return sequences


def blast_align(genome_sequence, target_sequence):
    """Alinea una secuencia de consulta contra un archivo FASTA usando BLASTn y devuelve un DataFrame."""

    # Ejecutar BLASTn con subprocess
    command_blast = [
        "blastn",
        "-query", genome_sequence,
        "-subject", target_sequence,
        "-outfmt", "6 sacc qstart qend"
    ]

    try:
        result = subprocess.run(command_blast, check=True, capture_output=True, text=True)

        # Revisar si hay salida vacía (sin alineaciones)
        output_lines = result.stdout.strip().split("\n")
        if not output_lines or output_lines[0] == "":
            print("No se encontraron alineaciones en BLAST.")
            return None

        # Convertir la salida en un DataFrame
        df = pd.DataFrame([line.split("\t") for line in output_lines], columns=["sacc", "qstart", "qend"])
        return df

    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar BLAST: {e}")
        return None


        # Convertir a DataFrame
        df = pd.DataFrame([line.split("\t") for line in output_lines], columns=["sacc", "qstart", "qend"])
        return df

    except subprocess.CalledProcessError as e:
        print(f"Error al ejecutar BLAST: {e}")
        return None


def gc_content(sequence_file):
    """
    Evalúa el contenido GC de una secuencia de nucleótidos.

    Parámetros:
        sequence_file (str): Puede ser una ruta de archivo FASTA o una cadena de texto con la secuencia.

    Retorna:
        float: Porcentaje de contenido GC.
    """
    total_length = 0
    gc_count = 0

    # Verificar si sequence_file es una cadena de texto que representa una secuencia
    if isinstance(sequence_file, str) and not sequence_file.endswith(('.fasta', '.fa', '.txt')):
        # Si es una cadena de texto, procesarla directamente
        sequence = sequence_file.upper()  # Convertir a mayúsculas para evitar problemas con minúsculas
        total_length = len(sequence)
        gc_count = sequence.count('G') + sequence.count('C')
    else:
        # Si es un archivo, procesarlo línea por línea
        try:
            with open(sequence_file) as file:
                for line in file:
                    if line.startswith('>'):
                        continue  # Ignorar la cabecera
                    line = line.strip().upper()  # Eliminar saltos de línea y convertir a mayúsculas
                    total_length += len(line)
                    gc_count += line.count('G') + line.count('C')
        except FileNotFoundError:
            return 0  # Si el archivo no existe, retornar 0

    if total_length == 0:
        return 0  # Evitar división por cero si la secuencia está vacía

    gc_percentage = round((gc_count / total_length) * 100, 2)
    return gc_percentage



def extract_range_fasta(secuence_file, start, end):
    """
    Extrae un rango específico de una secuencia en un archivo FASTA sin cargar toda la secuencia en memoria.

    Parámetros:
        file_path (str): Ruta al archivo FASTA.
        start (int): Posición inicial del rango (base 1).
        end (int): Posición final del rango (base 1).

    Retorna:
        str: Fragmento de la secuencia correspondiente al rango [start, end].
    """
    fragment = ""
    current_position = 0  # Contador de posición actual en la secuencia

    with open(secuence_file) as file:
        for line in file:
            if line.startswith('>'):
                continue  # Ignorar la primera linea

            line = line.strip()  # Eliminar saltos de línea
            line_length = len(line)

            # Calcular las posiciones relativas al fragmento que queremos extraer
            fragment_start = max(0, start - 1 - current_position)  # Ajustar a índice base 0
            fragment_end = max(0, end - current_position)

            # Si la línea actual contiene parte del fragmento que buscamos
            if fragment_start < line_length and fragment_end > 0:
                fragment += line[fragment_start:fragment_end]

            # Actualizar la posición actual
            current_position += line_length

            # Si hemos pasado el rango deseado, salir del bucle
            if current_position >= end:
                break

    return fragment


def find_pam_sites(dna_sequence, pam="NGG"):
    """ Encuentra todos los sitios PAM en la secuencia de ADN. """
    pam_regex = pam.replace("N", "[ATCG]")
    matches = [(m.start(), m.group()) for m in re.finditer(pam_regex, dna_sequence)]
    return matches




def design_sgRNAs_nn(gen_file,target_file, window_size=20):
    """ Genera candidatos a sgRNA en base a la región del target y sitios PAM. """

    target_ = blast_align(gen_file,target_file)
    qstart = target_.at[0, 'qstart']
    qend = target_.at[0, 'qend']

    qstart = int(qstart) - 20
    qend = int(qend) + 20

    target_region = extract_range_fasta(gen_file,qstart,qend)
    target_region_rvr_cmpl = complemento_inverso(target_region)

    pam_sites_rvr_cmpl = find_pam_sites(target_region_rvr_cmpl)


    pam_sites = find_pam_sites(target_region)
    df_sgRNA = pd.DataFrame(columns=["gRNA", "PAM", "GC_content", "position"])

    for site, pam_seq in pam_sites:
        if site >= window_size:
            candidate = target_region[site - window_size: site]
            gc_content_ = gc_content(candidate)
            eficiencia_ = round(predecir_eficiencia_nn(candidate),2)
            hebra = "+"
            #print(candidate)

            if 40 <= gc_content_ <= 80:
                new_row = pd.DataFrame([{
                    "gRNA": candidate,
                    "PAM": pam_seq,
                    "GC_content": gc_content_,
                    "position": site - window_size,
                    "hebra": hebra,
                    "Eficiencia": eficiencia_
                }])
                df_sgRNA = pd.concat([df_sgRNA, new_row], ignore_index=True)


    for site, pam_seq in pam_sites_rvr_cmpl:
        if site >= window_size:
            candidate = target_region_rvr_cmpl[site - window_size: site]
            gc_content_ = gc_content(candidate)
            eficiencia_ = round(predecir_eficiencia_nn(candidate),2)
            hebra = "-"
            #print(candidate)

            if 40 <= gc_content_ <= 80:
                new_row = pd.DataFrame([{
                    "gRNA": candidate,
                    "PAM": pam_seq,
                    "GC_content": gc_content_,
                    "position": site - window_size,
                    "hebra": hebra,
                    "Eficiencia": eficiencia_
                }])
                df_sgRNA = pd.concat([df_sgRNA, new_row], ignore_index=True)
    return df_sgRNA



def design_sgRNAs_xgb(gen_file,target_file, window_size=20):
    """ Genera candidatos a sgRNA en base a la región del target y sitios PAM. """

    target_ = blast_align(gen_file,target_file)
    qstart = target_.at[0, 'qstart']
    qend = target_.at[0, 'qend']

    qstart = int(qstart) - 20
    qend = int(qend) + 20

    target_region = extract_range_fasta(gen_file,qstart,qend)
    target_region_rvr_cmpl = complemento_inverso(target_region)

    pam_sites_rvr_cmpl = find_pam_sites(target_region_rvr_cmpl)


    pam_sites = find_pam_sites(target_region)
    df_sgRNA = pd.DataFrame(columns=["gRNA", "PAM", "GC_content", "position"])

    for site, pam_seq in pam_sites:
        if site >= window_size:
            candidate = target_region[site - window_size: site]
            gc_content_ = gc_content(candidate)
            eficiencia_ = round(predecir_eficiencia_xgb(candidate),2)
            hebra = "+"
            #print(candidate)

            if 40 <= gc_content_ <= 80:
                new_row = pd.DataFrame([{
                    "gRNA": candidate,
                    "PAM": pam_seq,
                    "GC_content": gc_content_,
                    "position": site - window_size,
                    "hebra": hebra,
                    "Eficiencia": eficiencia_
                }])
                df_sgRNA = pd.concat([df_sgRNA, new_row], ignore_index=True)


    for site, pam_seq in pam_sites_rvr_cmpl:
        if site >= window_size:
            candidate = target_region_rvr_cmpl[site - window_size: site]
            gc_content_ = gc_content(candidate)
            eficiencia_ = round(predecir_eficiencia_xgb(candidate),2)
            hebra = "-"
            #print(candidate)

            if 40 <= gc_content_ <= 80:
                new_row = pd.DataFrame([{
                    "gRNA": candidate,
                    "PAM": pam_seq,
                    "GC_content": gc_content_,
                    "position": site - window_size,
                    "hebra": hebra,
                    "Eficiencia": eficiencia_
                }])
                df_sgRNA = pd.concat([df_sgRNA, new_row], ignore_index=True)
    return df_sgRNA



def design_sgRNAs_rf(gen_file,target_file, window_size=20):
    """ Genera candidatos a sgRNA en base a la región del target y sitios PAM. """

    target_ = blast_align(gen_file,target_file)
    qstart = target_.at[0, 'qstart']
    qend = target_.at[0, 'qend']

    qstart = int(qstart) - 20
    qend = int(qend) + 20

    target_region = extract_range_fasta(gen_file,qstart,qend)
    target_region_rvr_cmpl = complemento_inverso(target_region)

    pam_sites_rvr_cmpl = find_pam_sites(target_region_rvr_cmpl)


    pam_sites = find_pam_sites(target_region)
    df_sgRNA = pd.DataFrame(columns=["gRNA", "PAM", "GC_content", "position"])

    for site, pam_seq in pam_sites:
        if site >= window_size:
            candidate = target_region[site - window_size: site]
            gc_content_ = gc_content(candidate)
            eficiencia_ = round(predecir_eficiencia(candidate),2)
            hebra = "+"
            #print(candidate)

            if 40 <= gc_content_ <= 80:
                new_row = pd.DataFrame([{
                    "gRNA": candidate,
                    "PAM": pam_seq,
                    "GC_content": gc_content_,
                    "position": site - window_size,
                    "hebra": hebra,
                    "Eficiencia": eficiencia_
                }])
                df_sgRNA = pd.concat([df_sgRNA, new_row], ignore_index=True)


    for site, pam_seq in pam_sites_rvr_cmpl:
        if site >= window_size:
            candidate = target_region_rvr_cmpl[site - window_size: site]
            gc_content_ = gc_content(candidate)
            eficiencia_ = round(predecir_eficiencia(candidate),2)
            hebra = "-"
            #print(candidate)

            if 40 <= gc_content_ <= 80:
                new_row = pd.DataFrame([{
                    "gRNA": candidate,
                    "PAM": pam_seq,
                    "GC_content": gc_content_,
                    "position": site - window_size,
                    "hebra": hebra,
                    "Eficiencia": eficiencia_
                }])
                df_sgRNA = pd.concat([df_sgRNA, new_row], ignore_index=True)
    return df_sgRNA
