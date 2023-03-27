# -*- coding: utf-8 -*-
"""
Module for de/encoding nucleotide sequences into tetrabinary strings

Author: ian.michael.bollinger@gmail.com    
"""
# Import Libraries and Set Globals
from PIL import Image
from sequence_diff import sequence_diff
import numpy as np
import bitarray
import os

global sequence_description
global nucleotide_sequence
global final_decoded_sequence
global decoded_description

def ttb_encode(nucleotide_sequence):
    """
    Encodes a given nucleotide sequence into a binary string using tetrabin encoding.

    Args:
    - nucleotide_sequence (str): A string representing the nucleotide sequence to be encoded.

    Returns:
    - final_encoded_string (str): A binary string representing the encoded nucleotide sequence.
    - nucleotide_type (str): A string representing the type of nucleotide used in the input sequence (DNA or RNA).

    Example usage:
    >>> ttb_encode("ATCG")
    ('100001100001001000010100001000', 'DNA')
    >>> ttb_encode("AUCG")
    ('100000000001001000010000001000', 'RNA')
    """

    # Define the tetrabin encoding scheme for nucleotides
    tetrabin_encoding_scheme = {'U': '0001', 'T': '0001', 'A': '0010', 'C': '0100', 'G': '1000', 
                                'R': '1010', 'Y': '0101', 'K': '1001', 'M': '0110',
                                'S': '1100', 'W': '0011', 'B': '1101', 'D': '1011',
                                'H': '0111', 'V': '1110', 'N': '1111', ' ': '0000'}
    
    # Define the spacer encoding scheme for nucleotide types
    spacer_encoding = {'DNA': '100001100001', 'RNA': '100000000001'}
    
    # Convert the nucleotide sequence to upper case and remove any new line characters
    nucleotide_sequence = nucleotide_sequence.upper()
    nucleotide_sequence = nucleotide_sequence.replace('\n','')
    
    # Determine the type of nucleotide (DNA or RNA) in the sequence
    if 'U' in nucleotide_sequence:
        nucleotide_type = 'RNA'
    else:
        nucleotide_type = 'DNA'
    
    # Encode the nucleotide sequence using tetrabin encoding
    try:
        encoded_sequence = ''.join(tetrabin_encoding_scheme[n] for n in nucleotide_sequence)
    except KeyError:
        print(f'Invalid nucleotide sequence: {nucleotide_sequence}')
        return None
    
    # Concatenate the encoded nucleotide sequence with the spacer for the corresponding nucleotide type
    encoded_spacer = spacer_encoding[nucleotide_type]
    final_encoded_string = encoded_spacer + encoded_sequence
    
    return(final_encoded_string, nucleotide_type)


def ttb_decode(final_encoded_string):
    """
    Decode a tetrabin-encoded sequence and return the decoded description, nucleotide type and sequence.

    Args:
    final_encoded_string (str): Tetrabin-encoded string to be decoded.

    Returns:
    tuple: A tuple containing the decoded description, nucleotide type and sequence.

    Raises:
    ValueError: If no spacer is found in the encoded binary sequence or if the nucleotide type is invalid.
    """

    # Tetrabin decoding scheme
    tetrabin_decoding_scheme = {'0001': 'T', '0010': 'A', '0100': 'C', '1000': 'G', 
                                '1010': 'R', '0101': 'Y', '1001': 'K', '0110': 'M', 
                                '1100': 'S', '0011': 'W', '1101': 'B', '1011': 'D',
                                '0111': 'H', '1110': 'V', '1111': 'N', '0000': ' '}

    # Spacer decoding scheme
    spacer_decoding = {'100001100001': 'DNA', '100000000001': 'RNA'}

    # Default nucleotide type
    nucleotide_type = None

    # Split the encoded sequence into description and sequence using the spacer
    for spacer in spacer_decoding.keys():
        if spacer in final_encoded_string:
            nucleotide_type = spacer_decoding[spacer]
            parts = final_encoded_string.split(spacer, 1)
            separated_description = parts[0]
            separated_sequence = parts[1]
            break

    # Raise an error if no spacer is found
    if separated_sequence is None:
        raise ValueError("No spacer found in encoded binary sequence.")

    # Decode the description and sequence
    decoded_description = "".join([chr(int(separated_description[i:i+8], 2)) for i in range(0, len(separated_description), 8)])
    decoded_sequence = ''.join(tetrabin_decoding_scheme[separated_sequence[i:i+4]] for i in range(0, len(separated_sequence), 4))

    # Replace 'T' with 'U' in the decoded sequence if it's RNA
    if nucleotide_type == 'RNA':
        final_decoded_sequence = decoded_sequence.replace('T','U')
    elif nucleotide_type == 'DNA':
        final_decoded_sequence = decoded_sequence
    else:
        raise ValueError("Invalid nucleotide type.")

    return (decoded_description, nucleotide_type, final_decoded_sequence)


def utf8_bin_encode(input_string):
    """
    Encodes the given input string using the UTF-8 encoding and returns the binary representation.

    Args:
    input_string (str): The input string to be encoded.

    Returns:
    str: The binary representation of the UTF-8 encoded input string.
    """

    # Encode input string using UTF-8
    utf8_bytes = input_string.encode('utf-8')

    # Convert UTF-8 encoded bytes to binary string
    output_binary = ''.join(format(byte, '08b') for byte in utf8_bytes)

    return output_binary


def utf8_bin_decode(input_binary):
    """
    Decodes the given binary string to its original string representation using the UTF-8 encoding.

    Args:
    input_binary (str): The binary string to be decoded.

    Returns:
    str: The original string representation of the decoded binary string.
    """

    # Check if input is already a binary string
    if all(c in ('0', '1') for c in input_binary):
        # If input is binary, skip conversion steps and directly decode the input
        byte_array = bitarray.bitarray(input_binary).tobytes()
        output_string = byte_array.decode('utf-8')
    else:
        # Split binary string into groups of 8 bits
        byte_list = [input_binary[i:i+8] for i in range(0, len(input_binary), 8)]
        
        # Convert each group of 8 bits to its decimal value
        decimal_list = [int(byte, 2) for byte in byte_list]
        
        # Convert the decimal values to their corresponding bytes
        byte_array = bytes(decimal_list)
        
        # Decode the bytes using the UTF-8 encoding to obtain the original input string
        output_string = byte_array.decode('utf-8')

    return output_string


def ascii_bin_encode(input_string: str) -> bitarray.bitarray:
    """
    Encodes a string of ASCII characters into a bitarray of binary-encoded ASCII characters.
    
    Args:
    - input_string: A string of ASCII characters to be encoded.
    
    Returns:
    - A bitarray of binary-encoded ASCII characters.
    """
    # Convert each character in the input string to its corresponding 8-bit ASCII binary code
    encoded_ascii_bin = ''.join(format(ord(char), '08b') for char in input_string)
    
    # Convert the binary-encoded ASCII characters into a bitarray
    ascii_bitarray = bitarray.bitarray(encoded_ascii_bin)
    
    return ascii_bitarray


def ascii_bin_decode(input_bitarray: bitarray.bitarray) -> str:
    """
    Decodes a bitarray of binary-encoded ASCII characters into a string of ASCII characters.
    
    Args:
    - input_bitarray: A bitarray of binary-encoded ASCII characters to be decoded.
    
    Returns:
    - A string of decoded ASCII characters.
    """
    # Convert the bitarray of binary-encoded ASCII characters to a binary string
    ASCII_Sequence_bin = input_bitarray.to01()

    # Split the binary string into groups of 8 bits
    byte_list = [ASCII_Sequence_bin[i:i+8] for i in range(0, len(ASCII_Sequence_bin), 8)]
    
    # Convert each group of 8 bits back to its decimal value
    decimal_list = [int(byte, 2) for byte in byte_list]
    
    # Convert each decimal value to its corresponding ASCII character
    char_list = [chr(decimal) for decimal in decimal_list]
    
    # Concatenate all the ASCII characters to form the original input string
    output_string = ''.join(char_list)
    
    return output_string

    
def string_to_list(s, group_size):
    """
    Convert a binary string to a list of lists, where each inner list
    contains 'group_size' number of integers from the string. If the length
    of the string is not a multiple of group_size, the last inner list will
    contain fewer elements.

    Args:
    - s (str): A binary string.
    - group_size (int): The number of bits to include in each inner list.

    Returns:
    - A list of lists containing integers.
    """
    return [list(map(int, s[i:i+group_size])) for i in range(0, len(s), group_size)]


def list_to_image4(l):
    """
    Convert a list of lists containing binary values to a PIL Image object,
    where each 0 in the list is represented by a black pixel and each 1 is
    represented by a white pixel. The image is duplicated 4 times in order to
    produce an RGBA image.

    Args:
    - l (list): A list of lists containing binary values.

    Returns:
    - A PIL Image object.
    """
    # Determine the dimensions of the image
    height = len(l)
    width = max(map(len, l))

    # Create a new image with a white background
    background_color = (255, 255, 255)
    new_image = Image.new('1', (width, height), background_color)

    # Add each pixel to the new image
    for i, row in enumerate(l):
        for j, value in enumerate(row):
            # Set the pixel color based on the value
            pixel_color = (0 if value else 255)

            # Set the pixel color for each of the four channels
            rgba = (pixel_color, pixel_color, pixel_color, pixel_color)
            new_image.putpixel((j, i), rgba)

    # Duplicate the image 4 times to produce an RGBA image
    rgba_image = new_image.convert('RGBA')
    r, g, b, a = rgba_image.split()
    rgba_image = Image.merge('RGBA', (r, g, b, a))

    return rgba_image


def fasta_encode(fasta_file_path):
    """
    Encode the contents of a FASTA file as a binary string, where the
    description line is UTF-8 encoded and concatenated with the nucleotide
    sequence, which is then compressed using the TTBinary algorithm.

    Args:
    - fasta_file_path (str): The path to the input FASTA file.

    Returns:
    - A tuple containing the encoded binary string and the original nucleotide
    sequence from the input file.
    """
    with open(fasta_file_path, 'r') as f:
        fasta_data = f.read()

    # Split the FASTA file into its components
    fasta_parts = fasta_data.split('\n', 1)
    sequence_description = fasta_parts[0].strip('>')
    nucleotide_sequence = fasta_parts[1]

    # Encode the description line as UTF-8
    description_binary = utf8_bin_encode(sequence_description)

    # Compress the nucleotide sequence using TTBinary
    final_encoded_string, nucleotide_type = ttb_encode(nucleotide_sequence)

    # Concatenate the description and compressed sequence
    final_encoded_string = description_binary + final_encoded_string

    return final_encoded_string, nucleotide_sequence


def png_encode(fasta_file_path):
    """
    Encodes a FASTA file into a PNG file using TTB encoding with RGBA pixel encoding.

    Args:
        fasta_file_path (str): The path to the input FASTA file.

    Returns:
        tuple: A tuple containing the path to the output PNG file and the nucleotide sequence from the input FASTA file.
    """
    # Create the path to the output PNG file
    fasta_name = fasta_file_path.split('/')[-1].replace('.','-')
    png_file_path = fasta_file_path.replace('.fasta', '.png').replace('FASTA FILES',f'OUTPUT FILES/{fasta_name}/') #png_file_path = fasta_file_path.replace('.fasta', '.png')
    
    # Encode the FASTA file using TTB encoding
    final_encoded_string, nucleotide_sequence = fasta_encode(fasta_file_path)
    
    # Encode the TTB-encoded string into an RGBA PNG image
    ttb_rgba_encoding(final_encoded_string, png_file_path)
    
    # Return the path to the output PNG file and the nucleotide sequence from the input FASTA file
    return png_file_path, nucleotide_sequence


def ttb_rgba_encoding(final_encoded_string, png_file_path):
    """
    Encodes a final RGBA string into a PNG image using the Text to Binary (TTB) encoding method.
    
    Args:
    - final_encoded_string (str): the final RGBA string to be encoded into a PNG image.
    - png_file_path (str): the path where the resulting PNG image should be saved.
    
    Raises:
    - ValueError: if the length of the input string is not a multiple of 4.
    
    Returns:
    - None
    """
    # Make sure the length of the input string is a multiple of 4
    if len(final_encoded_string) % 4 != 0:
        raise ValueError('FINAL RGBA STRING NOT DIVISIBLE BY 4')
    
    # Split the string into groups of four characters
    rgba_values = [final_encoded_string[i:i+4] for i in range(0, len(final_encoded_string), 4)]
    
    # Create a new image with a black background
    width, height = len(rgba_values), 1
    background_color = (0, 0, 0)
    new_image = Image.new('RGBA', (width, height), background_color)
    
    # Add each pixel to the new image
    for x, rgba in enumerate(rgba_values):
        # Convert the RGBA string into a list of integers
        rgba_list = [int(x) for x in list(rgba)]
        
        # Unpack the RGBA values from the list
        r, g, b, a = rgba_list
    
        # Add the pixel to the new image
        new_image.putpixel((x, 0), (r, g, b, a))
    
    # Save the image as a binary PNG file
    new_image.save(png_file_path, format='PNG')


def ttb_rgba_decoding(png_file_path):
    """
    Decodes a 1px-PNG image into a flattened binary string using the RGBA values of the pixel.

    Args:
    png_file_path (str): The file path to the PNG image file.

    Returns:
    final_encoded_string (str): The flattened binary string containing the decoded information from the PNG file.

    Raises:
    IOError: If the PNG file cannot be opened or decoded.
    """

    # Load the 1px-PNG image using Pillow
    try:
        loaded_img = Image.open(png_file_path)
    except IOError:
        raise IOError('Unable to open or decode PNG file: ' + png_file_path)

    # Extract the RGBA values of the pixel in the image
    img_data = np.asarray(loaded_img)
    flat_array = img_data.flatten()
    str_list = [str(item) for item in flat_array]

    # Combine the RGBA values into a flattened binary string
    delimiter = ''
    final_encoded_string = delimiter.join(str_list)

    return(final_encoded_string)


def png_decode(png_file_path):
    """
    Decode a PNG file that has been encoded using ttb_rgba_encoding.

    Args:
        png_file_path (str): Path to the PNG file to decode.

    Returns:
        A tuple containing the decoded description, nucleotide type, and final decoded sequence.
    """
    # Declare global variables to hold decoded data
    global final_decoded_sequence
    global decoded_description
    
    # Extract full encoded binary string from the PNG file
    final_encoded_string = ttb_rgba_decoding(png_file_path)
    
    # Decode the binary string to get the description, nucleotide type, and final decoded sequence
    decoded_description, nucleotide_type, final_decoded_sequence = ttb_decode(final_encoded_string)
    print(f'{nucleotide_type} SUCCESSFULLY DECODED PNG INTO DESCRIPTION, TYPE, & SEQUENCE')    
    
    return(decoded_description, nucleotide_type, final_decoded_sequence)


def reconstruct_fasta(decoded_description, nucleotide_type, final_decoded_sequence):
    """
    Reconstruct a FASTA-formatted string from the decoded description and sequence.

    Args:
        decoded_description (str): The description extracted from the encoded file.
        nucleotide_type (str): The type of nucleotide (e.g., DNA or RNA).
        final_decoded_sequence (str): The final decoded sequence.

    Returns:
        A string containing the reconstructed FASTA-formatted data.
    """
    reconstructed_fasta = f'>{decoded_description} {nucleotide_type}\n{final_decoded_sequence}'
    return(reconstructed_fasta)

 
# Debug & Example Area
if __name__ == '__main__':    
    working_directory = 'C:/Users/theda/OneDrive/Documents/Python/TetraBin-Nucleotide-Encoding'
    # Set fasta file path
    fasta_file_path = f'{working_directory}/Resources/FASTA FILES/AY281023.1.fasta'
    fasta_name = fasta_file_path.split('/')[-1]
    fasta_path = fasta_file_path.strip(fasta_name)

    temp_fasta_name = fasta_name.replace('.','-')
    fasta_folder_path = fasta_path.replace('FASTA FILES','OUTPUT FILES') + temp_fasta_name
    if not os.path.exists(fasta_folder_path):
        os.makedirs(fasta_folder_path)
       
    # Encode FASTA to PNG
    final_encoded_string, nucleotide_sequence = fasta_encode(fasta_file_path)
        
    png_file_path, nucleotide_sequence = png_encode(fasta_file_path)
    
    # Clean and prepare nucleotide sequence for comparison
    nucleotide_sequence = nucleotide_sequence.upper()
    nucleotide_sequence = nucleotide_sequence.replace('\n','')
    
    # Decode PNG to FASTA
    decoded_description, nucleotide_type, final_decoded_sequence = png_decode(png_file_path)
    
    # Compare original and decoded nucleotide sequence
    sequence_diff(nucleotide_sequence, final_decoded_sequence, 'Nucleotide Sequence')
    
    # Encode and Decode using BitArrays
    with open(fasta_file_path, 'r') as f:
        raw_fasta=f.read()
    
    # Encode ASCII
    encoded_ascii_bitarray = ascii_bin_encode(raw_fasta)
    decoded_fasta_bitarray = ascii_bin_decode(encoded_ascii_bitarray)
    
    # Compare original and decoded ASCII encoded FASTA
    sequence_diff(raw_fasta, decoded_fasta_bitarray, 'ASCII FASTA Data')
    
    # Encode UTF-8
    encoded_utf8_bitarray = bitarray.bitarray(utf8_bin_encode(raw_fasta))        
    decoded_utf8_bitarray = utf8_bin_decode(encoded_utf8_bitarray.to01())
    
    # Compare original and decoded UTF-8 encoded FASTA
    sequence_diff(raw_fasta, decoded_utf8_bitarray, 'UTF-8 FASTA Data')
