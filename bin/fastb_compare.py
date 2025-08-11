#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  2 08:06:47 2025

Binary Nucleotide Encoding (FASTB -> 'Fast-Binary')

@author: EYE
"""
from fastb import fasta_conversion
import os, subprocess, logging, zlib, bz2, time, argparse, shutil
from bitarray import bitarray
from rich.logging import RichHandler
import pandas as pd
from Bio import SeqIO
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = 'browser'

# Logging configuration
log_file = "compressor.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        RichHandler(rich_tracebacks=True, markup=True),
        logging.FileHandler(log_file, mode='w')
    ]
)
logger = logging.getLogger("compressor")


def get_file_size(file_path):
    """
    Get the size of a file in bytes.

    Parameters
    ----------
    file_path : str
        Path to the file.

    Returns
    -------
    int
        Size of the file in bytes.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    return os.path.getsize(file_path)


def load_fasta_to_df(fasta_path):
    """
    Load sequences from a FASTA file into a pandas DataFrame.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: 'id', 'description', 'sequence'.

    Raises
    ------
    FileNotFoundError
        If the specified file does not exist.
    """
    from Bio import SeqIO

    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"File not found: {fasta_path}")

    records = list(SeqIO.parse(fasta_path, "fasta"))
    data = [{
        "id": rec.id,
        "description": rec.description,
        "sequence": str(rec.seq)
    } for rec in records]

    return pd.DataFrame(data)


def generate_output_path(input_path, method):
    """
    Generate the appropriate output filename for compression.

    Parameters
    ----------
    input_path : str
        Path to the original input file.
    method : str
        Compression method ('pigz', 'zlib', or 'bz2').

    Returns
    -------
    str
        Output path with method-specific suffix.
    
    Raises
    ------
    ValueError
        If the method is unsupported.
    """
    base, _ = os.path.splitext(input_path)
    if method == 'pigz':
        return f"{base}_pigz.gz"
    elif method == 'zlib':
        return f"{base}_zlib.gz"
    elif method == 'bz2':
        return f"{base}.bz2"
    else:
        raise ValueError(f"Unsupported method: {method}")


def pigz_file(input_file_path, threads=4, decompress=False, preserve=False):
    """
    Compress or decompress a file using pigz (parallel gzip).

    Parameters
    ----------
    input_file_path : str
        Path to the file to compress or decompress.
    threads : int, optional
        Number of threads to use during compression (default is 4).
    decompress : bool, optional
        If True, decompress the file instead of compressing.
    preserve : bool, optional
        If True, the original file is preserved (copied instead of renamed). Default is True.

    Returns
    -------
    str
        Path to the output file.

    Raises
    ------
    ValueError
        If the input file extension is not '.gz' when decompressing.
    subprocess.CalledProcessError
        If the pigz command fails.
    """
    try:
        start = time.perf_counter()

        if decompress:
            if not input_file_path.endswith(".gz"):
                raise ValueError("Pigz decompression requires a '.gz' file.")
            output_file_path = input_file_path[:-3]
            logger.info(f"[cyan]Decompressing (pigz):[/] {input_file_path} → {output_file_path}")
            subprocess.run(["pigz", "-d", "-f", input_file_path], check=True)

        else:
            output_file_path = generate_output_path(input_file_path, 'pigz')
            temp_input = output_file_path.replace('_pigz.gz', '')

            # ✅ Preserve the original file
            if preserve:
                shutil.copy(input_file_path, temp_input)
            else:
                os.rename(input_file_path, temp_input)

            logger.info(f"[green]Compressing (pigz):[/] {input_file_path} → {output_file_path}")
            subprocess.run(["pigz", "-p", str(threads), "-f", temp_input], check=True)
            os.rename(temp_input + ".gz", output_file_path)

        duration = time.perf_counter() - start
        logger.info(f"[bold green]✔ Done:[/] {output_file_path} ({get_file_size(output_file_path)} bytes) in {duration:.2f}s")
        return output_file_path

    except Exception as e:
        logger.exception("Pigz compression/decompression failed.")
        raise e


def zlib_file(input_file_path, decompress=False, preserve=False):
    """
    Compress or decompress a file using zlib, saving output as a gzip-compatible .gz file.

    Parameters
    ----------
    input_file_path : str
        Path to the file to compress or decompress.
    decompress : bool, optional
        If True, decompress the file (default is False).

    Returns
    -------
    str
        Path to the output file.

    Raises
    ------
    ValueError
        If decompressing a file not ending in '.gz'.
    zlib.error
        If compression or decompression fails.
    """
    try:
        start = time.perf_counter()

        if decompress:
            if not input_file_path.endswith(".gz"):
                raise ValueError("Decompression requires a '.gz' file.")
            output_file_path = input_file_path[:-3]
            logger.info(f"[cyan]Decompressing (zlib):[/] {input_file_path} → {output_file_path}")
            with open(input_file_path, 'rb') as f_in, open(output_file_path, 'wb') as f_out:
                d = zlib.decompressobj(wbits=zlib.MAX_WBITS | 16)
                while chunk := f_in.read(1024):
                    f_out.write(d.decompress(chunk))
                f_out.write(d.flush())
        else:
            output_file_path = generate_output_path(input_file_path, 'zlib')            
            logger.info(f"[green]Compressing (zlib):[/] {input_file_path} → {output_file_path}")
            with open(input_file_path, 'rb') as f_in, open(output_file_path, 'wb') as f_out:
                c = zlib.compressobj(wbits=zlib.MAX_WBITS | 16)
                while chunk := f_in.read(1024):
                    f_out.write(c.compress(chunk))
                f_out.write(c.flush())

        duration = time.perf_counter() - start
        logger.info(f"[bold green]✔ Done:[/] {output_file_path} ({get_file_size(output_file_path)} bytes) in {duration:.2f}s")
        return output_file_path

    except Exception as e:
        logger.exception("Zlib compression/decompression failed.")
        raise e


def bz2_file(input_file_path, decompress=False, preserve=False):
    """
    Compress or decompress a file using BZ2 format.

    Parameters
    ----------
    input_file_path : str
        Path to the file to compress or decompress.
    decompress : bool, optional
        If True, decompress the file (default is False).

    Returns
    -------
    str
        Path to the output file.

    Raises
    ------
    ValueError
        If decompressing a file not ending in '.bz2'.
    OSError
        If reading/writing fails.
    """
    try:
        start = time.perf_counter()

        if decompress:
            if not input_file_path.endswith('.bz2'):
                raise ValueError("Decompression requires a '.bz2' file.")
            output_file_path = input_file_path[:-4]
            logger.info(f"[cyan]Decompressing (bz2):[/] {input_file_path} → {output_file_path}")
            with bz2.BZ2File(input_file_path, 'rb') as f_in, open(output_file_path, 'wb') as f_out:
                f_out.write(f_in.read())
        else:
            output_file_path = generate_output_path(input_file_path, 'bz2')                
            logger.info(f"[green]Compressing (bz2):[/] {input_file_path} → {output_file_path}")
            with open(input_file_path, 'rb') as f_in, bz2.BZ2File(output_file_path, 'wb') as f_out:
                f_out.write(f_in.read())

        duration = time.perf_counter() - start
        logger.info(f"[bold green]✔ Done:[/] {output_file_path} ({get_file_size(output_file_path)} bytes) in {duration:.2f}s")
        return output_file_path

    except Exception as e:
        logger.exception("BZ2 compression/decompression failed.")
        raise e


def plot_file_size_comparison(file_sizes, ncf_dict, fasta_name="FASTA File", color_scheme=None):
    """
    Plot a bar graph comparing file sizes (in bytes) and Nucleotide Compression Factor (NCF).

    Parameters
    ----------
    file_sizes : dict
        Dictionary of file sizes in bytes with labels as keys.
    ncf_dict : dict
        Dictionary of NCF values with same keys as file_sizes.
    fasta_name : str, optional
        Title for the plot, usually the FASTA filename.
    color_scheme : list of str, optional
        Optional list of colors for bars.

    Returns
    -------
    plotly.graph_objs._figure.Figure
        A Plotly bar chart.
    """
    # Compile into DataFrame
    df = pd.DataFrame({
        "Method": list(file_sizes.keys()),
        "File Size (B)": list(file_sizes.values()),
        "NCF (bp/byte)": [ncf_dict[k] for k in file_sizes.keys()]
    })

    # Create label string for each bar
    df["Label"] = df.apply(lambda row: f"{row['File Size (B)']:,} B<br>NCF: {row['NCF (bp/byte)']:.2f}", axis=1)

    # Use color scheme if provided
    colors = (
        color_scheme
        if color_scheme and len(color_scheme) == len(df)
        else ["darkred", "#FF7518", "orange", "#FFAE42", "green", "teal", "#0D98BA","turquoise"]
    )

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=df["Method"],
        y=df["File Size (B)"],
        text=df["Label"],
        textposition='inside',
        marker=dict(color=colors)
    ))

    fig.update_layout(
        title=dict(
            text=f"{fasta_name}<br>File Size & Nucleotide Compression Factor (NCF) Comparison",
            y=0.95,
            x=0.5,
            xanchor='center',
            yanchor='top',
            font=dict(size=30)
        ),
        yaxis_title=dict(
            text='File Size (bytes)',
            font=dict(size=20)
        ),
        xaxis_tickfont=dict(size=15),
        showlegend=False
    )

    fig.update_traces(textfont_size=13)
    return fig


if __name__ == "__main__":
    print("""
                          \033[0m,;\033[94m%%%%%%%%%,   \033[0m,;\033[94m%%%%%%%,   \033[0m,;\033[94m%%%%%%%, \033[0m,;\033[94m%%%%%%%%%%%, \033[0m,;\033[94m%%%%%%%%%%,
  \033[96m███████████\033[93m███████████  \033[0m;;\033[94m%%%%%%%%%% \033[0m;;\033[94m%%%%%%%%%%% \033[0m;;\033[94m%%%%%%%%% \033[0m;;\033[94m%%%%%%%%%%%% \033[0m;;\033[94m%%%%%%%%%%%
  \033[96m███████████\033[93m███████████  \033[0m;%\033[94m@@@\033[0m''''''' \033[0m;%\033[94m@@@\033[0m''';%\033[94m@@@ \033[0m;%\033[94m@@@\033[0m'''''   '''';%\033[94m@@@\033[0m'''' ;%\033[94m@@@\033[0m'''''\033[94m###
  \033[96m███████████\033[93m███████████  \033[0m;%\033[94m@@@        \033[0m;%\033[94m@@@   \033[0m;%\033[94m@@@ \033[0m;%\033[94m@@@            \033[0m;%\033[94m@@@     \033[0m;%\033[94m@@@    \033[0m;%\033[94m###
  \033[96m███████████\033[93m███████████  \033[0m;@\033[94m#########  \033[0m;@\033[94m########### \033[0m;@\033[94m#####,         \033[0m;@\033[94m###     \033[0m;@\033[94m###   \033[0m;@\033[94m###
  \033[96m███████████\033[93m███████████  \033[0m;@\033[94m#########  \033[0m;@\033[94m###########  \033[0m;@\033[94m#######,,     \033[0m;@\033[94m###     \033[0m;@\033[94m##########
  \033[92m███████████\033[91m███████████  \033[0m;@\033[94m###\033[0m;;;;;'  \033[0m;@\033[94m###\033[0m;;;;\033[94m####   \033[0m';;@\033[94m######;    \033[0m;@\033[94m###     \033[0m;@\033[94m###########
  \033[92m███████████\033[91m███████████  \033[0m;@\033[94m###        \033[0m;@\033[94m###    ####      \033[0m';;@\033[94m####    \033[0m;@\033[94m###     \033[0m;@\033[94m###\033[0m;;;;;;\033[94m###
  \033[92m███████████\033[91m███████████  \033[0m;@\033[94m###        \033[0m;@\033[94m###    ####       \033[0m';@\033[94m###'    \033[0m;@\033[94m###     \033[0m;@\033[94m###     \033[0m;@\033[94m###
  \033[92m███████████\033[91m███████████  \033[0m;#\033[94m888        \033[0m;#\033[94m888    8888 \033[0m;#\033[94m8888888888     \033[0m;#\033[94m888     \033[0m;#\033[94m888888888888
  \033[92m███████████\033[91m███████████  \033[0m;#\033[94m$$$        \033[0m;#\033[94m$$$    $$$$  \033[0m`#\033[94m$$$$$$$$      \033[0m;#\033[94m$$$     \033[0m;#\033[94m$$$$$$$$$$#
                              \033[94m╔════════════════════════════════════════════════════════════╗
                              \033[94m║\033[0m          Fast-Binary Nucleotide Encoder & Decoder          \033[94m║
                              \033[94m╚════════════════════════════════════════════════════════════╝
  
                         Curated & Maintained by Ian M Bollinger               
                              (\033[94mian.bollinger@entheome.org)\033[0m
  
                                   fastb_compare.py
                                      version  2
    """)

    parser = argparse.ArgumentParser(
        description="Fast-Binary Nucleotide Encoding & Decoding Comparison"
    )
    
    # Establish Default Values for fallback
    default_input_fasta = "/run/media/EYE/Main/TESTING_SPACE/nextflow_test/Psilocybe_subtropicalis/Psilocybe_subtropicalis-Entheome/Psilocybe_subtropicalis-Entheome_final_curated.fasta" # "/var/home/EYE/Desktop/Python Programs/FASTB/OR140556.fasta"
    default_cpu_threads = 12
    default_compare = True
    
    # Parse Input Arguments    
    parser.add_argument(
        "-i", "--input",
        type=str,
        default=default_input_fasta,
        help="Path to the input FASTA file"
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=default_cpu_threads,
        help=f"Number of CPU threads to use for pigz compression (default: {default_cpu_threads})"
    )
    parser.add_argument(
        "-c", "--compare",
        action="store_true",
        default=default_compare,
        help=f"Generate compression comparison graph (default: {default_compare})"
    )
    args = parser.parse_args()
    input_fasta = args.input
    cpu_threads = args.threads
    compare_bool = args.compare

    # Provide Logger input information
    logger.info(f"[bold yellow]Input File:[/] {input_fasta}")
    logger.info(f"[bold yellow]Threads:[/] {cpu_threads}")
    logger.info(f"[bold yellow]Compare Graph:[/] {compare_bool}")

    # Load FASTA Data
    try:
        fasta_df = load_fasta_to_df(input_fasta)
        total_bp = sum(len(rec.seq) for rec in SeqIO.parse(input_fasta, "fasta"))
        logger.info(f"[bold green]✔ Loaded FASTA:[/] {len(fasta_df)} sequences from '{input_fasta}'")
        logger.info(f"[bold cyan]Total base pairs across all sequences:[/] {total_bp:,} bp")
    except Exception as e:
        logger.error(f"[red]Failed to load FASTA file:[/] {e}")
        raise

    # CONVERT FASTA DATA INTO FASTB and encode into .fastb data file
    converted_fastb = fasta_conversion(input_fasta)

    # COMPARISON GRAPH
    if compare_bool:

        # TODO: CONVERT input_fasta INTO ASCII binary data file

        # TODO: CONVERT input_fasta INTO UTF8 binary data file       

        # Typical Compression
        fasta_pigz = pigz_file(input_fasta, threads=cpu_threads, preserve=compare_bool)
        fasta_zlib = zlib_file(input_fasta, preserve=compare_bool)
        fasta_bz2 = bz2_file(input_fasta, preserve=compare_bool)
        fastb_pigz = pigz_file(converted_fastb, threads=cpu_threads, preserve=compare_bool)
        fastb_zlib = zlib_file(converted_fastb, preserve=compare_bool)
        fastb_bz2 = bz2_file(converted_fastb, preserve=compare_bool)

        file_sizes = {
            "fasta": get_file_size(input_fasta),
            # TODO: add ascii_bin
            # TODO: add utf8_bin
            "fasta_pigz": get_file_size(fasta_pigz),
            "fasta_zlib": get_file_size(fasta_zlib),
            "fasta_bz2": get_file_size(fasta_bz2),
            "fastb": get_file_size(converted_fastb),
            "fastb_pigz": get_file_size(fastb_pigz),
            "fastb_zlib": get_file_size(fastb_zlib),
            "fastb_bz2": get_file_size(fastb_bz2),
        }
        
        # Compute Nucleotide Compression Factor (NCF)
        ncf_dict = {k: round(total_bp / v, 2) for k, v in file_sizes.items()}
        logger.info("[bold magenta]Nucleotide Compression Factor (bp/byte):[/]")
        for method, ratio in ncf_dict.items():
            logger.info(f"  {method}: {ratio} bp/byte")        
        
        # Generate a Comparison Plot
        comparison_fig = plot_file_size_comparison(file_sizes, ncf_dict, fasta_name=os.path.basename(input_fasta))
        output_svg = input_fasta.replace(".fasta","_fastb_comp.svg")
        comparison_fig.write_image(output_svg)
        comparison_fig.show()
        
