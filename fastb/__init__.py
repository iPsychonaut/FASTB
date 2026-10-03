"""FASTB: 2-bit binary sequence format (file format 3.1)."""
from fastb.core import (
    Record,
    write_file,
    read_file,
    read_index,
    read_record,
    iter_raw,
    encode_sequence,
    decode_sequence,
    MAGIC_LINE,
)
from fastb.io import iter_records, write_records, iter_fasta
from fastb.alphabet import (
    AlphabetError,
    ProteinDetectedError,
    detect_alphabet,
)

__version__ = "3.1.0"
__all__ = [
    "Record", "write_file", "read_file", "read_index", "read_record",
    "MAGIC_LINE", "iter_records", "write_records", "iter_fasta", "iter_raw", "encode_sequence", "decode_sequence",
    "AlphabetError", "ProteinDetectedError", "detect_alphabet",
    "__version__",
]
