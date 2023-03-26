# -*- coding: utf-8 -*-
"""
Module for comparing two sequences of the same length.

Author: ian.michael.bollinger@gmail.com

Functions:
    sequence_diff(seq1, seq2, seq_type)
        Compares two sequences of the same length and prints the differing indices.

Arguments:
    seq1: str
        The first sequence to be compared.
    seq2: str
        The second sequence to be compared.
    seq_type: str
        A string describing the type of sequence being compared.

Raises:
    ValueError:
        If the length of seq1 and seq2 are not the same.

Returns:
    sequence_validation: True/False

Example Usage:
    desc_seq1 = '10000.1 Organism description'
    desc_seq2 = '10000.1 Organism description'
    sequence_diff(desc_seq1, desc_seq2, 'Description')

    nuc_seq1 = 'TACG'
    nuc_seq2 = 'TACG'
    sequence_diff(nuc_seq1, nuc_seq2, 'Nucleotide Sequence')

    bin_seq1 = '0101001001001000'
    bin_seq2 = '0001001001001000'
    sequence_diff(bin_seq1, bin_seq2, 'Binary Sequence')
    
"""

# Compares two sequences of the same length and prints the differing indices.
def sequence_diff(seq1, seq2, seq_type):
    # Establish Default Validation State
    sequence_validation = False
    
    # Check that the sequences are of the same length
    if len(seq1) != len(seq2):
        raise ValueError(f"{seq_type}s are not the same length")
    
    # Initialize list to store differing indices
    differing_indices = []

    # Iterate over the sequences and compare the characters at each index
    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if a != b:
            differing_indices.append(i)
    
    # Change Validation State if needed Return and Print the result
    if not differing_indices:
        print(f"{seq_type}s are identical, QC PASS")
        sequence_validation = True
    else:
        print(f"{seq_type}s differ at index/indices: {differing_indices}")
    return(sequence_validation)
        
if __name__ == '__main__':
    # Example string containing the Sequence's Description
    desc_seq1 = '10000.1 Organism description'
    desc_seq2 = '10000.1 Organism description'
    sequence_diff(desc_seq1, desc_seq2, 'Description')
    
    # Example string containing the Nucleotide Sequence
    nuc_seq1 = 'TACG'
    nuc_seq2 = 'TACC'
    sequence_diff(nuc_seq1, nuc_seq2, 'Nucleotide Sequence')

    # Example string containing the Binary Sequence
    bin_seq1 = '001001001001000'
    bin_seq2 = '0001001001001000'
    sequence_diff(bin_seq1, bin_seq2, 'Binary Sequnce')