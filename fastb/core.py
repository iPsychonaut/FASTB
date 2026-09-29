"""FASTB v3 reference implementation.

A FASTB file stores nucleotide sequences packed at 2 bits per base (4 bits
for records with IUPAC codes other than N). Everything about the sequence
(magic line, record headers, index) is plain UTF-8 text. See FASTB-SPEC.md
at the repo root for the normative layout.

Layout at a glance:

    ##FASTB 3.0                              <- magic line, byte 0
    # optional file comments
    >chr1  LEN=12  ENC=2  ALPHA=D  CRC=...  BYTES=3  NRUNS=...  MASK=...
    <3 payload bytes>
    \\n
    >chr2  ...
    chr1  85   12                            <- index entries (plain text)
    chr2  157  12
    ##INDEX  ENTRIES=2  OFFSETS_START=...
    ##END    FILE_CRC=...
"""

from __future__ import annotations
import zlib
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import numpy as np

from fastb.alphabet import require_nucleotide

MAGIC_LINE = b"##FASTB 3.0\n"

# 2-bit: A=00 C=01 G=10 T/U=11. N packs as 00 and is restored from NRUNS.
_PACK2 = np.full(256, 255, dtype=np.uint8)
for _ch, _code in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("U", 3), ("N", 0)):
    _PACK2[ord(_ch)] = _code

# 4-bit IUPAC one-hot: bit 3=A, bit 2=C, bit 1=G, bit 0=T/U.
_IUPAC4 = {
    "A": 0b1000, "C": 0b0100, "G": 0b0010, "T": 0b0001, "U": 0b0001,
    "W": 0b1001, "S": 0b0110, "M": 0b1100, "K": 0b0011,
    "R": 0b1010, "Y": 0b0101,
    "B": 0b0111, "D": 0b1011, "H": 0b1101, "V": 0b1110,
    "N": 0b1111,
    "-": 0b0000, ".": 0b0000,
}
_PACK4 = np.full(256, 255, dtype=np.uint8)
for _ch, _code in _IUPAC4.items():
    _PACK4[ord(_ch)] = _code

# Uppercase table: a-z -> A-Z, everything else unchanged.
_UPPER = np.arange(256, dtype=np.uint8)
_UPPER[ord("a"):ord("z") + 1] -= 32


def _unpack_tables(alpha: str):
    """Byte -> ASCII lookup tables for one alphabet. Shape (256, 4) and (256, 2)."""
    t = "T" if alpha == "D" else "U"
    sym2 = np.frombuffer(("ACG" + t).encode(), dtype=np.uint8)
    b = np.arange(256, dtype=np.uint8)
    lut2 = np.stack([sym2[(b >> s) & 3] for s in (6, 4, 2, 0)], axis=1)

    sym4 = np.full(16, ord("-"), dtype=np.uint8)
    for ch, code in _IUPAC4.items():
        if ch not in ("U", "."):
            sym4[code] = ord(ch)
    if alpha == "R":
        sym4[0b0001] = ord("U")
    lut4 = np.stack([sym4[b >> 4], sym4[b & 0xF]], axis=1)
    return lut2, lut4


_UNPACK = {"D": _unpack_tables("D"), "R": _unpack_tables("R")}


class Record:
    """A single nucleotide record in a FASTB v3 file."""

    __slots__ = ("name", "sequence", "alpha", "comment")

    def __init__(self, name: str, sequence: str, alpha: str = "D",
                 comment: Optional[str] = None):
        if "\t" in name or " " in name or "\n" in name or name.startswith(">"):
            raise ValueError("name must not contain whitespace, newline, or '>'")
        if alpha not in ("D", "R"):
            raise ValueError(
                f"alpha must be 'D' (DNA) or 'R' (RNA); got {alpha!r}. "
                f"('P' is reserved for protein in a future release.)"
            )
        self.name = name
        self.sequence = sequence
        self.alpha = alpha
        self.comment = comment

    def __repr__(self):
        return f"Record(name={self.name!r}, len={len(self.sequence)}, alpha={self.alpha})"


# ---------------------------------------------------------------------------
# Run lists: (start, length) pairs, hex, comma-separated in the header
# ---------------------------------------------------------------------------

def _runs(flags: np.ndarray) -> List[Tuple[int, int]]:
    """(start, length) of each contiguous True run in a bool array."""
    if not flags.any():
        return []
    # Positions where the flag changes value. Kept as bool ops: an int8
    # diff with Python-int endpoints upcasts to int64 (8x the memory).
    change = np.flatnonzero(flags[1:] != flags[:-1]) + 1
    starts = change[flags[change]]
    ends = change[~flags[change]]
    if flags[0]:
        starts = np.concatenate(([0], starts))
    if flags[-1]:
        ends = np.concatenate((ends, [len(flags)]))
    return list(zip(starts.tolist(), (ends - starts).tolist()))


def _encode_runs(runs: List[Tuple[int, int]]) -> str:
    return ",".join(f"{s:x}:{l:x}" for s, l in runs)


def _decode_runs(s: str) -> List[Tuple[int, int]]:
    if not s:
        return []
    return [(int(a, 16), int(b, 16)) for a, b in (p.split(":") for p in s.split(","))]


# ---------------------------------------------------------------------------
# Codec
# ---------------------------------------------------------------------------

def encode_sequence(seq: str, alpha: str) -> Tuple[int, bytes, str, str]:
    """Pack one sequence. Returns (enc, payload, nruns_field, mask_field)."""
    raw = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    up = _UPPER[raw]
    mask = _encode_runs(_runs(raw != up))

    codes = _PACK2[up]
    if not (codes == 255).any():
        nruns = _encode_runs(_runs(up == ord("N")))
        n = len(codes)
        padded = np.zeros(-(-n // 4) * 4, dtype=np.uint8)
        padded[:n] = codes
        q = padded.reshape(-1, 4)
        packed = (q[:, 0] << 6) | (q[:, 1] << 4) | (q[:, 2] << 2) | q[:, 3]
        return 2, packed.astype(np.uint8).tobytes(), nruns, mask

    codes = _PACK4[up]
    if (codes == 255).any():
        bad = chr(int(up[codes == 255][0]))
        raise ValueError(f"illegal symbol {bad!r} for ALPHA={alpha}")
    n = len(codes)
    padded = np.zeros(-(-n // 2) * 2, dtype=np.uint8)
    padded[:n] = codes
    q = padded.reshape(-1, 2)
    packed = (q[:, 0] << 4) | q[:, 1]
    return 4, packed.astype(np.uint8).tobytes(), "", mask


def decode_sequence(payload: bytes, length: int, enc: int, alpha: str,
                    nruns: List[Tuple[int, int]], mask: List[Tuple[int, int]]) -> str:
    """Unpack one payload back to a sequence string."""
    lut2, lut4 = _UNPACK[alpha]
    b = np.frombuffer(payload, dtype=np.uint8)
    if enc == 2:
        out = lut2[b].reshape(-1)[:length]
    elif enc == 4:
        out = lut4[b].reshape(-1)[:length]
    else:
        raise ValueError(f"Unknown ENC={enc}")
    out = np.ascontiguousarray(out)
    for s, l in nruns:
        out[s:s + l] = ord("N")
    for s, l in mask:
        out[s:s + l] |= 0x20
    return out.tobytes().decode("ascii")


# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------

def write_file(records: Iterable[Record], out, file_comment: Optional[str] = None,
               force_nucleotide: bool = False) -> None:
    """Write records to `out` (binary file-like with .write() and .tell()).

    Records are written as they arrive; only the index (a few bytes per
    record) is held in memory until the end.
    """
    out.write(MAGIC_LINE)
    if file_comment:
        for line in file_comment.splitlines():
            out.write(b"# " + line.encode("utf-8") + b"\n")

    index_entries: List[Tuple[str, int, int]] = []

    for rec in records:
        require_nucleotide(rec.sequence, rec.alpha, rec.name, force_nucleotide)
        enc, payload, nruns, mask = encode_sequence(rec.sequence, rec.alpha)
        crc = zlib.crc32(payload) & 0xFFFFFFFF

        if rec.comment:
            comment_safe = rec.comment.replace("\n", " ").replace("\t", " ")
            out.write(b"# " + comment_safe.encode("utf-8") + b"\n")

        index_entries.append((rec.name, out.tell(), len(rec.sequence)))

        header = (
            f">{rec.name}\tLEN={len(rec.sequence)}\tENC={enc}\tALPHA={rec.alpha}"
            f"\tCRC={crc:08x}\tBYTES={len(payload)}"
        )
        if nruns:
            header += f"\tNRUNS={nruns}"
        if mask:
            header += f"\tMASK={mask}"
        out.write(header.encode("utf-8") + b"\n")
        out.write(payload)
        out.write(b"\n")

    index_start = out.tell()
    for name, off, length in index_entries:
        out.write(f"{name}\t{off}\t{length}\n".encode("utf-8"))
    out.write(
        f"##INDEX\tENTRIES={len(index_entries)}\tOFFSETS_START={index_start}\n".encode("utf-8")
    )
    # Streaming writers cannot know the file CRC up front; 00000000 means "not set".
    out.write(b"##END\tFILE_CRC=00000000\n")


# ---------------------------------------------------------------------------
# Read
# ---------------------------------------------------------------------------

def _parse_header_line(line: bytes) -> Dict[str, str]:
    if not line.startswith(b">") or not line.endswith(b"\n"):
        raise ValueError(f"Malformed record header: {line!r}")
    parts = line[1:-1].decode("utf-8").split("\t")
    fields = {"NAME": parts[0]}
    for p in parts[1:]:
        if "=" not in p:
            raise ValueError(f"Malformed header field (no '='): {p!r}")
        k, v = p.split("=", 1)
        fields[k] = v
    return fields


def _record_from_fields(fields: Dict[str, str], payload: bytes,
                        comment: Optional[str] = None) -> Record:
    name = fields["NAME"]
    alpha = fields.get("ALPHA")
    if alpha is None:
        raise ValueError(f"Record {name!r}: missing ALPHA field in header.")
    if alpha not in ("D", "R"):
        raise ValueError(
            f"Record {name!r}: ALPHA={alpha} is not supported by this fastb "
            f"version. Known values: D (DNA), R (RNA). "
            f"(P is reserved for protein in a future release.)"
        )
    expected_crc = int(fields["CRC"], 16)
    actual_crc = zlib.crc32(payload) & 0xFFFFFFFF
    if actual_crc != expected_crc:
        raise ValueError(
            f"CRC mismatch for {name!r}: header says {expected_crc:08x}, "
            f"payload computes {actual_crc:08x}"
        )
    seq = decode_sequence(
        payload, int(fields["LEN"]), int(fields["ENC"]), alpha,
        _decode_runs(fields.get("NRUNS", "")), _decode_runs(fields.get("MASK", "")),
    )
    return Record(name=name, sequence=seq, alpha=alpha, comment=comment)


def iter_raw(src, read_payload: bool = True) -> Iterator[Tuple[Dict[str, str], bytes, Optional[str]]]:
    """Walk records without decoding. Yields (header_fields, payload, comment).

    With read_payload=False the payload is skipped with seek() and b'' is
    yielded in its place; the CRC is not checked.
    """
    magic = src.readline()
    if magic != MAGIC_LINE:
        raise ValueError(f"Not a FASTB v3 file. Expected magic {MAGIC_LINE!r}, got {magic!r}")

    pending_comment: Optional[str] = None
    while True:
        line = src.readline()
        if not line:
            return
        if line.startswith(b"#"):
            if line.startswith(b"##INDEX") or line.startswith(b"##END"):
                return
            pending_comment = line[1:].strip().decode("utf-8", errors="replace")
            continue
        if not line.startswith(b">"):
            if line.strip() == b"":
                continue
            if b"\t" in line:
                return  # first index entry: past the last record
            raise ValueError(f"Unexpected line outside record: {line!r}")

        fields = _parse_header_line(line)
        byte_len = int(fields["BYTES"])
        if read_payload:
            payload = src.read(byte_len)
            if len(payload) != byte_len:
                raise ValueError(
                    f"Truncated payload for {fields['NAME']!r}: "
                    f"expected {byte_len} bytes, got {len(payload)}"
                )
        else:
            payload = b""
            src.seek(byte_len, 1)
        term = src.read(1)
        if term != b"\n":
            raise ValueError(
                f"Missing record terminator after {fields['NAME']!r} "
                f"(got {term!r}, expected b'\\n')"
            )
        yield fields, payload, pending_comment
        pending_comment = None


def read_file(src) -> Iterator[Record]:
    """Iterate decoded records from `src` (binary file-like). Streams one record at a time."""
    for fields, payload, comment in iter_raw(src):
        yield _record_from_fields(fields, payload, comment)


# ---------------------------------------------------------------------------
# Random access via footer index
# ---------------------------------------------------------------------------

def read_index(path: str) -> Dict[str, Tuple[int, int]]:
    """Read the footer index. Returns {name: (header_offset, length_in_bases)} in file order."""
    with open(path, "rb") as f:
        f.seek(0, 2)
        file_size = f.tell()
        back = min(file_size, 4096)
        f.seek(file_size - back)
        tail = f.read(back)

        idx_marker = tail.rfind(b"\n##INDEX\t")
        if idx_marker < 0:
            raise ValueError("No ##INDEX trailer found")
        f.seek(file_size - back + idx_marker + 1)
        idx_line = f.readline().decode("utf-8")
        parts = dict(p.split("=", 1) for p in idx_line.strip().split("\t")[1:])
        offsets_start = int(parts["OFFSETS_START"])
        entries = int(parts["ENTRIES"])

        f.seek(offsets_start)
        result = {}
        for _ in range(entries):
            name, off, length = f.readline().decode("utf-8").strip().split("\t")
            result[name] = (int(off), int(length))
        return result


def read_record(path: str, key) -> Record:
    """Read one record by name (str) or by 0-based position (int), using the footer index."""
    idx = read_index(path)
    if isinstance(key, int):
        if not 0 <= key < len(idx):
            raise KeyError(f"Record index {key} out of range (file has {len(idx)} records)")
        header_offset = list(idx.values())[key][0]
    else:
        if key not in idx:
            raise KeyError(f"Record {key!r} not in file")
        header_offset = idx[key][0]
    with open(path, "rb") as f:
        f.seek(header_offset)
        fields = _parse_header_line(f.readline())
        payload = f.read(int(fields["BYTES"]))
        return _record_from_fields(fields, payload)
