"""FASTB v3: binary sequence format."""
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
from fastb.alphabet import (
    AlphabetError,
    ProteinDetectedError,
    detect_alphabet,
)

__version__ = "3.0.0"
__all__ = [
    "Record", "write_file", "read_file", "read_index", "read_record",
    "MAGIC_LINE", "iter_raw", "encode_sequence", "decode_sequence",
    "AlphabetError", "ProteinDetectedError", "detect_alphabet",
    "__version__",
]
