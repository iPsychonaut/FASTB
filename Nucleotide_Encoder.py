"""
Module for de/encoding nucleotide sequences into tetrabinary strings

Author: ian.michael.bollinger@gmail.com
    
"""
###############################################################################
# IMPORTS & SETTINGS AREA
import bz2
import bitarray
import configparser
import gzip
import io
import os
import struct
import sys
import zlib
import numpy as np
import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go
from Bio import SeqIO
from PIL import Image
from plotly.offline import plot

# Use Python's built-in configparser library to parse the variables in the config.txt file
config = configparser.ConfigParser()
config.read('C:/Users/theda/OneDrive/Documents/Python/ttb_config.txt') # REPLACE WITH YOUR CONFIG LOCATION

modules_path = config.get('DEFAULT', 'modules_path')
fasta_path = config.get('DEFAULT', 'fasta_path')
output_path = config.get('DEFAULT', 'output_path')

# Append the directory containing Custom Modules to sys.path
sys.path.append(modules_path)

# Declare Default Plotly Renderer
pio.renderers.default='svg'

# Import Custom Modules
from sequence_diff import sequence_diff
from nucleotide_tetrabin import (
    ascii_bin_encode, 
    ascii_bin_decode, 
    utf8_bin_encode, 
    utf8_bin_decode,
    fasta_encode,
    png_encode, 
    png_decode, 
    reconstruct_fasta
)


###############################################################################
# FUCTION AREA
def get_file_size(file_path):
    """
    This function returns the size of a file in bytes.

    Args:
    file_path (str): The path of the file whose size is to be determined.

    Returns:
    int: The size of the file in bytes.
    """
    return os.path.getsize(file_path)


def compress_gzip(input_file_path, output_file_path):
    """
    This function compresses a file using GZIP compression.

    Args:
    input_file_path (str): The path of the input file to be compressed.
    
    output_file_path (str): The path of the compressed output file.

    Returns:
    str: The path of the compressed file.
    """
    with open(input_file_path, 'rb') as f_in:
        with gzip.open(output_file_path, 'wb') as f_out:
            f_out.write(f_in.read())
    return output_file_path


def compress_zlib(input_file_path, output_file_path):
    """
    This function compresses a file using ZLIB compression.

    Args:
    input_file (str): The path of the input file to be compressed.
    
    output_file_path (str): The path of the compressed output file.

    Returns:
    str: The path of the compressed file.
    """
    with open(input_file_path, 'rb') as f_in:
        with open(output_file_path, 'wb') as f_out:
            compressor = zlib.compressobj()
            while True:
                chunk = f_in.read(1024)
                if not chunk:
                    break
                compressed_chunk = compressor.compress(chunk)
                f_out.write(compressed_chunk)
            remaining_data = compressor.flush()
            f_out.write(remaining_data)
    return output_file_path


def compress_bz2(input_file_path, output_file_path):
    """
    This function compresses a file using BZ2 compression.

    Args:
    input_file (str): The path of the input file to be compressed.
    
    output_file_path (str): The path of the compressed output file.

    Returns:
    str: The path of the compressed file.
    """
    with open(input_file_path, 'rb') as f_in:
        with bz2.BZ2File(output_file_path, 'wb') as f_out:
            f_out.write(f_in.read())
    return output_file_path


def fasta_check(fasta_file_path,decoded_description,final_decoded_sequence):
    """
    Check if the reconstructed FASTA-formatted string matches the original FASTA file.

    Args:
        fasta_file_path (str): Path to the original FASTA file.
        decoded_description (str): The decoded description.
        final_decoded_sequence (str): The final decoded sequence.

    Returns:
        None.
    """
    fasta_df, input_fasta_data = fasta_to_df_data(fasta_file_path)

    reconstructed_fasta = f'>{decoded_description}\n{final_decoded_sequence}'
        
    # Check if the reconstructed FASTA-formatted string matches the original FASTA data
    sequence_diff(input_fasta_data, reconstructed_fasta, 'FASTA Data')


def fasta_to_df_data(fasta_file_path):
    # Read the FASTA file using Biopython's SeqIO.parse() function
    records = list(SeqIO.parse(fasta_file_path, "fasta"))
    
    # Convert the records to a pandas dataframe
    fasta_list = []
    for record in records:
        fasta_dict = {
            "Accession #": record.id,
            "Description": record.description.strip(record.id),
            "Nucleotide": str(record.seq)
        }
        fasta_list.append(fasta_dict)
        
    fasta_df = pd.DataFrame(fasta_list, index=range(len(fasta_list)), columns=["Accession #", "Description", "Nucleotide"])
    
    # Extract the first record from the dataframe and format it as a FASTA file string
    nucleotide_data = fasta_df.loc[0, 'Nucleotide']
    description_data = fasta_df.loc[0, 'Description']
    if 'U' in nucleotide_data.upper():
        nucleotide_type = 'RNA'
    else:
        if 'RNA' in description_data:
            nucleotide_type = 'RNA'
    input_fasta_data = f">{fasta_df.loc[0, 'Accession #']}{description_data}\n{nucleotide_data}"
    return(fasta_df, input_fasta_data)


###############################################################################
# INPUT AREA
# Define file paths
working_dir = fasta_path.strip('\Resources\FASTA FILES')
os.chdir(working_dir)

# fasta_file_path = f'{fasta_path}/AY281023.1.fasta'

fasta_file_path = f'{fasta_path}/NC_000016.10-89913418-89920973.fasta'

# fasta_file_path = f'{fasta_path}/CM039002.1.fasta'

fasta_name = fasta_file_path.split('/')[-1]
temp_fasta_name = fasta_name.replace('.','-')
fasta_folder_path = fasta_path.replace('\FASTA FILES','\OUTPUT FILES\\') + temp_fasta_name
if not os.path.exists(fasta_folder_path):
    os.makedirs(fasta_folder_path)
ascii_bin = fasta_folder_path + '//' + fasta_name.replace('.fasta', '_ascii.bin')
utf8_bin = fasta_folder_path + '//' + fasta_name.replace('.fasta', '_utf8.bin')
fasta_gzip_path = fasta_folder_path + '//' + fasta_name.replace('.fasta', '.fasta_gzip.gz')
fasta_zlib_path = fasta_folder_path + '//' + fasta_name.replace('.fasta', '.fasta_zlib.gz')
fasta_bz2_path = fasta_folder_path + '//' + fasta_name.replace('.fasta', '.fasta.bz2')
fasta_df, input_fasta_data = fasta_to_df_data(fasta_file_path)

###############################################################################
# ENCODING AREA
import array
import h5py

# FASTA file to Binary Tetrad 1-Pixel PNG
tetrabin_png, nucleotide_sequence = png_encode(fasta_file_path)
    
# FASTA String to ASCII binary
encoded_ascii_bitarray = ascii_bin_encode(input_fasta_data)

# Write bit sequence to bitarray file
with open(ascii_bin, 'wb') as f:
    encoded_ascii_bitarray.tofile(f)

# FASTA String to UTF-8 binary
encoded_utf8_bitarray = bitarray.bitarray(utf8_bin_encode(input_fasta_data))

# Write bit sequence to bitarray file
with open(utf8_bin, 'wb') as f:
    encoded_utf8_bitarray.tofile(f)


###############################################################################
# DECODE PNG
decoded_description, nucleotide_type, final_decoded_sequence = png_decode(tetrabin_png)


###############################################################################
# SEQUENCE QUALITY CONTROL CHECK
fasta_check(fasta_file_path,decoded_description,final_decoded_sequence)


###############################################################################
# COMPRESSION AREA           
fasta_gzip = compress_gzip(fasta_file_path, fasta_gzip_path)
fasta_zlib = compress_zlib(fasta_file_path, fasta_zlib_path)
fasta_bz2 = compress_bz2(fasta_file_path, fasta_bz2_path)


###############################################################################
# COMPILE FILE SIZES
file_sizes = {
    "fasta": get_file_size(fasta_file_path),
    "ascii_bin": get_file_size(ascii_bin),
    "utf8_bin": get_file_size(utf8_bin),
    "tetrabin_png": get_file_size(tetrabin_png),
    "fasta_gzip": get_file_size(fasta_gzip),
    "fasta_zlib": get_file_size(fasta_zlib),
    "fasta_bz2": get_file_size(fasta_bz2),
}


###############################################################################
# SIZE COMPARISON GRAPH AREA
x = list(file_sizes.keys())
y = list(file_sizes.values())
colors = ['darkred', 'red', 'orangered','green', 'lightskyblue', 'skyblue', 'lightblue']

# Create bar graph
fig = go.Figure()
fig.add_trace(go.Bar(x = x,
                      y = y,
                      text = y,
                      textposition = 'inside',
                      marker = dict(color=colors))) 

# update the layout of the graph
fig.update_layout(
    title={'text': f'{fasta_name}<br>File Size Comparison',
            'y': 0.95,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': 35}},
    yaxis_title={'text': 'File Size (bytes)',
                  'font': {'size': 25}},
    xaxis_tickfont={'size': 15},
    annotations=[dict(text='Original FASTA File',
                      x=0,
                      y=(file_sizes['fasta'] * 1.1),
                      font_size=15,
                      showarrow=False)] + 
                [dict(text=f"{round(((file_sizes['fasta']-file_sizes[x])/file_sizes['fasta'])*100,2)}% File<br>Size Reduction",
                      x=i,
                      y=(file_sizes[x] * 1.15),
                      font_size=15,
                      showarrow=False) for i, x in enumerate(x) if x != 'fasta'],
    showlegend=False)

# update the data labels
fig.update_traces(textposition='auto')

# Show bar graph
plot(fig)
