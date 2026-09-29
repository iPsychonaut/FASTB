"""Streaming helpers: FASTA in, FASTB out, and back. Nothing here buffers a whole file."""

from __future__ import annotations
import gzip
import os
from typing import Iterable, Iterator, Tuple

import numpy as np

from fastb.core import Record, read_file, write_file


def open_text_or_gz(path: str):
    """Open a FASTA file for binary reading. .gz is transparently decompressed."""
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def iter_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """Yield (name, sequence) per record. Name is the header up to the first space."""
    name = None
    parts = []
    with open_text_or_gz(path) as f:
        for line in f:
            line = line.rstrip(b"\r\n")
            if line.startswith(b">"):
                if name is not None:
                    seq = b"".join(parts).decode("ascii")
                    parts = []  # free the line list before the consumer encodes
                    yield name, seq
                header = line[1:].split()
                name = header[0].decode("utf-8") if header else ""
                parts = []
            elif line:
                parts.append(line)
    if name is not None:
        seq = b"".join(parts).decode("ascii")
        parts = []
        yield name, seq


def iter_records(path: str) -> Iterator[Tuple[str, str]]:
    """Yield (name, sequence) from a .fastb file, one record at a time."""
    with open(path, "rb") as f:
        for rec in read_file(f):
            yield rec.name, rec.sequence


def write_records(path: str, records: Iterable[Tuple[str, str]],
                  file_comment: str | None = None, force_nucleotide: bool = False) -> int:
    """Write (name, sequence) pairs to a .fastb file as they arrive. Returns record count.

    Writes to path + '.tmp' and renames on success, so a failed run leaves no
    partial .fastb behind. ALPHA is R if the record contains U, else D.
    """
    count = 0

    def gen():
        nonlocal count
        for name, seq in records:
            alpha = "R" if "U" in seq or "u" in seq else "D"
            count += 1
            yield Record(name, seq, alpha=alpha)

    tmp = path + ".tmp"
    try:
        with open(tmp, "wb") as f:
            write_file(gen(), f, file_comment=file_comment, force_nucleotide=force_nucleotide)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    return count


def wrap_bytes(seq: str, width: int) -> bytes:
    """Sequence as bytes with a newline every `width` characters and at the end."""
    raw = seq.encode("ascii")
    if width <= 0 or not raw:
        return raw + b"\n"
    arr = np.frombuffer(raw, dtype=np.uint8)
    n = len(arr)
    rows = -(-n // width)
    padded = np.full(rows * width, ord("\n"), dtype=np.uint8)
    padded[:n] = arr
    out = np.full((rows, width + 1), ord("\n"), dtype=np.uint8)
    out[:, :width] = padded.reshape(rows, width)
    data = out.tobytes()
    pad = rows * width - n  # unused cells in the last row
    if pad:
        data = data[:-(pad + 1)] + b"\n"
    return data
