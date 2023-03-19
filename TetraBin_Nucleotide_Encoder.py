# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 11:35:20 2023

@author: theda
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

pio.renderers.default='svg'

def compress_to_binary(nucleotide_sequence):
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
    
    # Convert the nucleotide sequence to binary
    binary_string = ''.join(encoding_scheme[n] for n in tqdm.tqdm(nucleotide_sequence, desc='Encoding', unit='nucleotides'))
    binary_sequence = binary_string.encode()
    return(binary_sequence)

def decompress_to_nucleotide(compressed_file):
    """
    Decodes a binary file to a nucleotide sequence using a custom encoding scheme.
    
    Parameters:
    - compressed_file: str
        The path of the binary file to decode.
    
    Returns:
    - decompressed_nucleotide: str
        The nucleotide sequence decoded from the binary file.
    """
    # Read the binary file
    with gzip.open(compressed_file, 'rb') as f:
        binary_sequence = f.read()
        
        # Define the reverse encoding scheme
        reverse_encoding_scheme = {'0000': 'T', '0001': 'A', '0010': 'C', '0011': 'G', 
                                  '0100': 'R', '0101': 'Y', '0110': 'K', '0111': 'M', 
                                  '1000': 'S', '1001': 'W', '1010': 'B', '1011': 'D',
                                  '1100': 'H', '1101': 'V', '1110': 'N'}
    
        # Detect the encoding of the binary data
        detected_encoding = chardet.detect(binary_sequence)['encoding']
        
        # If the detected encoding is None, raise an error
        if detected_encoding is None:
            raise ValueError('Could not detect encoding of binary data.')
        
        # Decode the binary data using the detected encoding
        binary_string = binary_sequence.decode(detected_encoding)
        
        # Convert the binary sequence to nucleotides
        nucleotide_sequence = []
        for i in tqdm.tqdm(range(0, len(binary_string), 4), desc='Decoding', unit='nucleotides'):
            nucleotide_sequence.append(reverse_encoding_scheme[binary_string[i:i+4]])
        
        decompressed_nucleotide = ''.join(nucleotide_sequence)
        return decompressed_nucleotide

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

###############################################################################
#
#
#
###############################################################################

# Open the FASTA file

fasta_file = 'C:/Users/theda/OneDrive/Documents/Python/AY281023.1.fasta'

# fasta_file = 'C:/Users/theda/OneDrive/Documents/Python/NC_000016.10-89913418-89920973.fasta'

# fasta_file = 'C:/Users/theda/OneDrive/Documents/Python/CM039002.1.fasta'

compressed_file = fasta_file.replace('.fasta','.bin')

# Use Biopython's SeqIO.parse() function to parse the FASTA file
records = list(SeqIO.parse(fasta_file, "fasta"))

# Create a pandas dataframe from the SeqIO records
df = pd.DataFrame(columns=["Accession #", "Description", "Nucleotide"])
data = []
for record in records:
    data.append({"Accession #": record.id, "Description": record.description, "Nucleotide": str(record.seq)})
df = pd.concat([df, pd.DataFrame(data)], ignore_index=True)

# Print the dataframe
print(df)

for row, sequence_data in df.iterrows():
    print(row)
    acc_numb = sequence_data['Accession #']
    description = sequence_data['Description']
    input_string = f'>{acc_numb} {description}'
    nucleotide_sequence = sequence_data['Nucleotide']
    binary_sequence = compress_to_binary(nucleotide_sequence)    
    df.loc[row, 'Binary'] = binary_sequence

    # Define the string and binary number sequence
    encoded_string = input_string.encode()
    string_len = len(encoded_string)
    binary_len = len(binary_sequence)
    
    # Use format codes to create a packed binary string
    packed_data = struct.pack(f"{string_len}s{binary_len}s", encoded_string, binary_sequence)
    
    # Compress the binary sequence
    compressed_sequence = gzip.compress(packed_data)
    
    # Check the compressed size
    compressed_size = len(compressed_sequence)
    print(f'Compressed size: {compressed_size}')
    
    # Decompress the binary sequence
    decompressed_sequence = gzip.decompress(compressed_sequence)
    
    sequence_diff(packed_data, decompressed_sequence)
    
    # Write the compressed data to a binary file
    with open(compressed_file, 'wb') as f:
        f.write(compressed_sequence)
        
# decompressed_nucleotide = decompress_to_nucleotide(compressed_file)

# sequence_diff(nucleotide_sequence, decompressed_nucleotide)


# get the file sizes in bytes
original_size = os.path.getsize(fasta_file)
compressed_size = os.path.getsize(compressed_file)

percent_difference_size = round((original_size-compressed_size)/original_size, 2) * 100

# create a bar graph
fig = go.Figure(data=[go.Bar(name='Original File',
                             x=['Original File'],
                             y=[original_size],
                             text=[original_size],
                             textfont=dict(size=26)),
                      go.Bar(name='Compressed File',
                             x=['Compressed File'],
                            y=[compressed_size],
                            text=[compressed_size],
                            textfont=dict(size=26))])

# update the layout of the graph
fig.update_layout(
    title={'text': 'File Sizes Comparison',
            'y': 0.95,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': 30}},
    yaxis_title={'text': 'File Size (bytes)',
                 'font': {'size': 20}},
    annotations=[
        dict(text=f'{percent_difference_size}% File Size Reduction',
             x=1,
             y=(compressed_size*1.5),
             font_size=32,
             showarrow=False)])

# update the data labels
fig.update_traces(textposition='auto')

# show the graph
plot(fig)