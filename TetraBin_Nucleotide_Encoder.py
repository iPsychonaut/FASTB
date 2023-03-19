# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 11:35:20 2023

@author: ian.michael.bollinger@gmail.com
"""

import pandas as pd
from Bio import SeqIO
import gzip
import struct
import chardet
import os
import time
import tqdm
import plotly 
import plotly.graph_objects as go
import plotly.io as pio
from plotly.offline import plot
import os

pio.renderers.default='svg'

def compress_to_binary(fasta_file_path):
    """
    Returns the binary sequence of a nucleotide sequence based on a custom encoding scheme.
    
    Parameters:
    - nucleotide_sequence: str
        The nucleotide sequence to be encoded.
        
    Returns:
    - binary_sequence: bytes
        The binary sequence of the nucleotide sequence.
    """
    # Define the custom encoding scheme
    encoding_scheme = {'U': '0000', 'T': '0000', 'A': '0001', 'C': '0010', 'G': '0011', 
                       'R': '0100', 'Y': '0101', 'K': '0110', 'M': '0111', 'S': '1000',
                       'W': '1001', 'B': '1010', 'D': '1011', 'H': '1100',
                       'V': '1101', 'N': '1110'}
    
    # Open the FASTA file
    fasta_name = fasta_file_path.split('/')[-1]

    compressed_file_path = fasta_file_path.replace('.fasta','.bin').replace('FASTA','BIN')

    # Use Biopython's SeqIO.parse() function to parse the FASTA file
    records = list(SeqIO.parse(fasta_file_path, "fasta"))

    # Create a pandas dataframe from the SeqIO records
    fasta_df = pd.DataFrame(columns=["Accession #", "Description", "Nucleotide"])
    data = []
    for record in records:
        fasta_dictionary = {"Accession #": record.id, "Description": record.description, "Nucleotide": str(record.seq)}
        data.append(fasta_dictionary)
    fasta_df = pd.concat([fasta_df, pd.DataFrame(data)], ignore_index=True)

    # Iterate over each the dataframe and generate a 
    for row, sequence_data in fasta_df.iterrows():
        acc_numb = sequence_data['Accession #']
        description = sequence_data['Description']
        input_string = f'>{acc_numb} {description} '
        nucleotide_sequence = sequence_data['Nucleotide']
        
        # Convert the nucleotide sequence to binary
        compressed_binary_string = ''.join(encoding_scheme[n] for n in tqdm.tqdm(nucleotide_sequence, desc='Encoding', unit=' nucleotides'))
        compressed_binary_sequence = compressed_binary_string.encode()
        fasta_df.loc[row, 'Binary'] = compressed_binary_sequence

        # Define the string and binary number sequence
        encoded_string = input_string.encode()
        string_len = len(encoded_string)
        binary_len = len(compressed_binary_sequence)
        
        # Use format codes to create a packed binary string
        packed_data = struct.pack(f"{string_len}s{binary_len}s", encoded_string, compressed_binary_sequence)
        
        # Compress the binary sequence
        compressed_sequence = gzip.compress(packed_data)
        
        # Decompress the binary sequence
        decompressed_sequence = gzip.decompress(compressed_sequence)
        
        # Check Differences
        sequence_diff(packed_data, decompressed_sequence)
        
        # Write the compressed data to a binary file
        with open(compressed_file_path, 'wb') as f:
            f.write(compressed_sequence)
            
    return(fasta_df, compressed_file_path)
    
def decompress_to_nucleotide(compressed_file_path):
    """
    Decodes a binary file to a nucleotide sequence using a custom encoding scheme.
    
    Parameters:
    - compressed_file_path: str
        The path of the binary file to decode.
    
    Returns:
    - decompressed_nucleotide: str
        The nucleotide sequence decoded from the binary file.
    """
    # Read the binary file
    with gzip.open(compressed_file_path, 'rb') as f:
        compressed_binary_sequence = str(f.read()).strip("b'")
        
    split = compressed_binary_sequence.split(' ')

    description = ''

    for i, item in enumerate(split):
        if i == 0:
            acc_numb = item.strip('>')
        elif i == len(split)-1:
            decompressed_binary_string = item
        else:
            description += f'{item} '

    # Define the reverse encoding scheme
    reverse_encoding_scheme = {'0000': 'T', '0001': 'A', '0010': 'C', '0011': 'G', 
                              '0100': 'R', '0101': 'Y', '0110': 'K', '0111': 'M', 
                              '1000': 'S', '1001': 'W', '1010': 'B', '1011': 'D',
                              '1100': 'H', '1101': 'V', '1110': 'N'}    


    # Convert the binary sequence to nucleotides
    nucleotide_sequence = []
    for i in tqdm.tqdm(range(0, len(decompressed_binary_string), 4), desc='Decoding', unit=' tetrads'):
        nucleotide_sequence.append(reverse_encoding_scheme[decompressed_binary_string[i:i+4]])

    decompressed_nucleotide = ''.join(nucleotide_sequence)

    decompressed_dictionary = {'Accession #': acc_numb,
                      'Description': description,
                      'Nucleotide': decompressed_nucleotide,
                      'Binary': decompressed_binary_string}
    decompressed_df = pd.DataFrame(decompressed_dictionary, index=[0])
    return(decompressed_df)

def sequence_diff(seq1, seq2):
    """
    Compares two sequences and returns the index/indices where they differ.

    Args:
        seq1 (str): First sequence to compare
        seq2 (str): Second sequence to compare

    Raises:
        ValueError: If the sequences are not of the same length.

    Returns:
        None: Prints the result to the console.

    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences are not the same length")
    
    # Initialize list to store differing indices
    differing_indices = []

    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if a != b:
            differing_indices.append(i)

    if not differing_indices:
        print("Sequences are identical")
    else:
        print(f"Sequences differ at index/indices: {differing_indices}")

def plot_file_size_comparison(fasta_file_path: str, compressed_file_path: str) -> None:
    """
    Compare the file size of the original and compressed files, and display the comparison as a bar graph.

    Args:
        fasta_file_path (str): Path to the original file.
        compressed_file_path (str): Path to the compressed file.

    Returns:
        None
    """
    # Fasta Name
    fasta_name = fasta_file_path.split('/')[-1]
    
    # get the file sizes in bytes
    original_size = os.path.getsize(fasta_file_path)
    compressed_size = os.path.getsize(compressed_file_path)

    percent_difference_size = round((original_size - compressed_size) / original_size, 2) * 100

    # create a bar graph
    fig = go.Figure(data=[
        go.Bar(name='Original File', x=['Original File'], y=[original_size], text=[original_size],
               textfont=dict(size=30)),
        go.Bar(name='Compressed File', x=['Compressed File'], y=[compressed_size], text=[compressed_size],
               textfont=dict(size=30))])

    # update the layout of the graph
    fig.update_layout(
        title={'text': f'{fasta_name}<br>File Size Comparison',
               'y': 0.95,
               'x': 0.5,
               'xanchor': 'center',
               'yanchor': 'top',
               'font': {'size': 35}},
        yaxis_title={'text': 'File Size (bytes)',
                     'font': {'size': 35}},
        xaxis_tickfont={'size': 35},
        annotations=[
            dict(text=f'{percent_difference_size}% File Size Reduction',
                 x=1,
                 y=(compressed_size * 1.5),
                 font_size=30,
                 showarrow=False)],
        showlegend=False)

    # update the data labels
    fig.update_traces(textposition='auto')

    # show the graph
    plot(fig)

###############################################################################
#
#
#
###############################################################################

# Get current working directory
working_directory = os.getcwd()

# fasta_file_path = f'{working_directory}/FASTA FILES/AY281023.1.fasta'

# fasta_file_path = f'{working_directory}/FASTA FILES/NC_000016.10-89913418-89920973.fasta'

fasta_file_path = f'{working_directory}/FASTA FILES/CM039002.1.fasta'

fasta_df, compressed_file_path = compress_to_binary(fasta_file_path)
        
decompressed_df = decompress_to_nucleotide(compressed_file_path)

sequence_diff(fasta_df['Nucleotide'][0], decompressed_df['Nucleotide'][0])

plot_file_size_comparison(fasta_file_path, compressed_file_path)