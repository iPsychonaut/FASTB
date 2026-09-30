"""Streaming helpers: FASTA in, FASTB out, and back. Nothing here buffers a whole file."""

from __future__ import annotations
import gzip
import math
import os
import zlib
from typing import Iterable, Iterator, List, Optional, Tuple

import numpy as np

from fastb.core import (CHUNK_BASES, Record, _decode_runs, _parse_header_line,
                        _run_bounds, decode_chunk, iter_raw, read_file, write_file)


def open_text_or_gz(path: str):
    """Open a FASTA file for binary reading. .gz is transparently decompressed."""
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def iter_fasta(path: str, as_bytes: bool = False) -> Iterator[Tuple[str, str]]:
    """Yield (name, sequence) per record. Name is the header up to the first space.

    Reads 8 MB blocks and strips newlines with numpy, so a 200 Mb record
    costs a few dozen 8 MB pieces instead of millions of line objects.
    as_bytes=True yields the sequence as bytes and saves one full copy when
    the record goes straight into the encoder.
    """
    name = None
    parts: List[bytes] = []
    buf = b""

    def flush():
        seq = b"".join(parts)
        parts.clear()
        return name, (seq if as_bytes else seq.decode("ascii"))

    with open_text_or_gz(path) as f:
        while True:
            block = f.read(1 << 23)
            if block:
                data = buf + block
                cut = data.rfind(b"\n")
                if cut < 0:
                    buf = data
                    continue
                data, buf = data[:cut + 1], data[cut + 1:]
            else:
                data, buf = buf, b""
                if not data:
                    break
            pos = 0
            while pos < len(data):
                if data[pos:pos + 1] == b">":
                    if name is not None:
                        yield flush()
                    eol = data.find(b"\n", pos)
                    if eol < 0:
                        eol = len(data)
                    header = data[pos + 1:eol].split()
                    name = header[0].decode("utf-8") if header else ""
                    pos = eol + 1
                else:
                    if name is None:
                        raise ValueError(f"{path}: sequence data before the first '>' header")
                    nxt = data.find(b"\n>", pos)
                    end = len(data) if nxt < 0 else nxt + 1
                    seg = np.frombuffer(data[pos:end], dtype=np.uint8)
                    parts.append(seg[(seg != 10) & (seg != 13)].tobytes())
                    pos = end
            if not block:
                break
    if name is not None:
        yield flush()


def iter_records(path: str) -> Iterator[Tuple[str, str]]:
    """Yield (name, sequence) from a .fastb file, one record at a time."""
    with open(path, "rb") as f:
        for rec in read_file(f):
            yield rec.name, rec.sequence


def write_records(path: str, records: Iterable[Tuple[str, str]],
                  file_comment: Optional[str] = None, force_nucleotide: bool = False) -> int:
    """Write (name, sequence) pairs to a .fastb file as they arrive. Returns record count.

    Writes to path + '.tmp' and renames on success, so a failed run leaves no
    partial .fastb behind. ALPHA is R if the record contains U, else D.
    """
    count = 0

    def gen():
        nonlocal count
        for name, seq in records:
            if isinstance(seq, str):
                alpha = "R" if "U" in seq or "u" in seq else "D"
            else:
                alpha = "R" if b"U" in seq or b"u" in seq else "D"
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


# ---------------------------------------------------------------------------
# FASTB -> FASTA, chunked, optionally across processes
# ---------------------------------------------------------------------------

def wrap_bytes(seq, width: int, final: bool = True) -> bytes:
    """ASCII bytes with a newline every `width` characters.

    `seq` is a str or a uint8 array. With final=True the output ends with a
    newline. With final=False the caller promises len(seq) is a multiple of
    width (or width is 0), so the output ends exactly at a row boundary.
    """
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8) if isinstance(seq, str) else seq
    n = len(arr)
    if width <= 0 or n == 0:
        return arr.tobytes() + (b"\n" if final else b"")
    rows = -(-n // width)
    padded = np.full(rows * width, ord("\n"), dtype=np.uint8)
    padded[:n] = arr
    out = np.full((rows, width + 1), ord("\n"), dtype=np.uint8)
    out[:, :width] = padded.reshape(rows, width)
    data = out.tobytes()
    pad = rows * width - n
    if pad:
        data = data[:-(pad + 1)] + b"\n"
    return data


def _chunk_bases(width: int, per_byte: int) -> int:
    """Largest chunk <= CHUNK_BASES that is a multiple of the wrap width and of bases-per-byte."""
    unit = per_byte if width <= 0 else per_byte * width // math.gcd(per_byte, width)
    return max(unit, CHUNK_BASES // unit * unit)


def _plan(fields: dict, payload_offset: int, width: int) -> List[tuple]:
    """Split one record into independent chunk tasks: (path-free) tuples for _decode_task."""
    enc = int(fields["ENC"])
    per_byte = 4 if enc == 2 else 2
    length = int(fields["LEN"])
    alpha = fields["ALPHA"]
    nruns = _run_bounds(_decode_runs(fields.get("NRUNS", "")))
    mask = _run_bounds(_decode_runs(fields.get("MASK", "")))
    step = _chunk_bases(width, per_byte)
    tasks = []
    for base in range(0, max(length, 1), step):
        count = min(step, length - base)
        nbytes = -(-count // per_byte)
        final = base + count >= length
        # Only the runs touching this chunk travel with the task (keeps pickles small).
        sub = []
        for starts, ends in (nruns, mask):
            lo = np.searchsorted(ends, base, side="right")
            hi = np.searchsorted(starts, base + count, side="left")
            sub.append((starts[lo:hi], ends[lo:hi]))
        tasks.append((payload_offset + base // per_byte, nbytes, enc, alpha, base, count,
                      sub[0], sub[1], width, final))
    return tasks


_FILES: dict = {}


def _decode_task(path: str, task: tuple) -> bytes:
    """Read one chunk of payload from `path` and return it as wrapped FASTA bytes."""
    offset, nbytes, enc, alpha, base, count, nruns, mask, width, final = task
    f = _FILES.get(path)
    if f is None:
        f = _FILES[path] = open(path, "rb")
    f.seek(offset)
    out = decode_chunk(f.read(nbytes), enc, alpha, base, count, nruns, mask)
    return wrap_bytes(out, width, final)


def _record_headers(path: str) -> Iterator[Tuple[dict, int]]:
    """Yield (header_fields, payload_offset) per record; verifies each payload CRC."""
    with open(path, "rb") as f:
        pos = 0
        for fields, _, _ in iter_raw(f, read_payload=False):
            # iter_raw left us just past the record terminator; recompute the payload start.
            end = f.tell()
            nbytes = int(fields["BYTES"])
            payload_offset = end - 1 - nbytes
            f.seek(payload_offset)
            crc = 0
            for i in range(0, nbytes, 1 << 24):
                crc = zlib.crc32(f.read(min(1 << 24, nbytes - i)), crc)
            if crc != int(fields["CRC"], 16):
                raise ValueError(
                    f"CRC mismatch for {fields['NAME']!r}: header says {fields['CRC']}, "
                    f"payload computes {crc:08x}")
            f.seek(end)
            yield fields, payload_offset


def stream_fasta(path: str, out, width: int = 80, threads: int = 1,
                 names: Optional[set] = None) -> int:
    """Write every record of a .fastb file to `out` as FASTA. Returns record count.

    Memory is bounded by CHUNK_BASES per worker. With threads > 1 the chunks
    are decoded by a process pool in file order while this process writes.
    """
    def tasks():
        for fields, payload_offset in _record_headers(path):
            if names is not None and fields["NAME"] not in names:
                continue
            yield fields["NAME"], _plan(fields, payload_offset, width)

    count = 0
    if threads <= 1:
        for name, plan in tasks():
            out.write(b">" + name.encode("utf-8") + b"\n")
            for t in plan:
                out.write(_decode_task(path, t))
            count += 1
        return count

    from multiprocessing import Pool
    from itertools import chain

    def flat():
        nonlocal count
        for name, plan in tasks():
            count += 1
            yield (name, None)
            for t in plan:
                yield (path, t)

    with Pool(threads) as pool:
        for item in pool.imap(_run_item, flat(), chunksize=1):
            out.write(item)
    return count


def _run_item(item) -> bytes:
    path_or_name, task = item
    if task is None:
        return b">" + path_or_name.encode("utf-8") + b"\n"
    return _decode_task(path_or_name, task)
