"""FASTB v3 reference implementation.

Nucleotide sequence container with compact bit-packed payload, human-readable
structural layer, explicit self-description, and O(1) random access by record
name. Clean break from v1 (sentinel bitstream) and v2 (TLV container).

File layout at a glance:

    ##FASTB 3.0                              <- magic line, byte 0
    # optional file comments
    >chr1  LEN=12  ENC=2  ALPHA=DNA  CRC=...  BYTES=3
    <3 payload bytes>
    \n
    >chr2  ...
    ...
    chr1  85   12                            <- index entries (plain text)
    chr2  157  12
    ##INDEX  ENTRIES=2  OFFSETS_START=...
    ##END    FILE_CRC=...

See FASTB_v3_spec.md for the normative specification.
"""

from __future__ import annotations
import io
import zlib
from typing import List, Tuple, Optional, Iterator, Dict

MAGIC_LINE = b"##FASTB 3.0\n"

# 2-bit: A=00 C=01 G=10 T/U=11
# A<->T and C<->G are bitwise complements, so reverse-complement is a byte NOT
# followed by a per-byte nibble swap.
_PACK2 = {"A": 0b00, "C": 0b01, "G": 0b10, "T": 0b11, "U": 0b11}
_UNPACK2_DNA = {0b00: "A", 0b01: "C", 0b10: "G", 0b11: "T"}
_UNPACK2_RNA = {0b00: "A", 0b01: "C", 0b10: "G", 0b11: "U"}

# 4-bit IUPAC one-hot: bit 3=A, bit 2=C, bit 1=G, bit 0=T/U.
# Degenerates are bitwise OR of the bases they represent.
_PACK4 = {
    "A": 0b1000, "C": 0b0100, "G": 0b0010, "T": 0b0001, "U": 0b0001,
    "W": 0b1001, "S": 0b0110, "M": 0b1100, "K": 0b0011,
    "R": 0b1010, "Y": 0b0101,
    "B": 0b0111, "D": 0b1011, "H": 0b1101, "V": 0b1110,
    "N": 0b1111,
    "-": 0b0000, ".": 0b0000,
}
_UNPACK4_DNA = {
    0b1000: "A", 0b0100: "C", 0b0010: "G", 0b0001: "T",
    0b1001: "W", 0b0110: "S", 0b1100: "M", 0b0011: "K",
    0b1010: "R", 0b0101: "Y",
    0b0111: "B", 0b1011: "D", 0b1101: "H", 0b1110: "V",
    0b1111: "N",
    0b0000: "-",
}
_UNPACK4_RNA = dict(_UNPACK4_DNA)
_UNPACK4_RNA[0b0001] = "U"

_LEGAL_DNA = set("ACGTWSMKRYBDHVN-.")
_LEGAL_RNA = set("ACGUWSMKRYBDHVN-.")


class Record:
    """A single nucleotide record in a FASTB v3 file."""

    __slots__ = ("name", "sequence", "alpha", "mask", "comment")

    def __init__(
        self,
        name: str,
        sequence: str,
        alpha: str = "DNA",
        mask: Optional[List[Tuple[int, int]]] = None,
        comment: Optional[str] = None,
    ):
        if "\t" in name or " " in name or "\n" in name or name.startswith(">"):
            raise ValueError("name must not contain whitespace, newline, or '>'")
        if alpha not in ("DNA", "RNA"):
            raise ValueError(
                f"alpha must be 'DNA' or 'RNA'; got {alpha!r}. "
                f"('AA' is reserved for a future release.)"
            )
        self.name = name
        self.sequence = sequence
        self.alpha = alpha
        self.mask = mask or []
        self.comment = comment

    def __repr__(self):
        return f"Record(name={self.name!r}, len={len(self.sequence)}, alpha={self.alpha})"


# ---------------------------------------------------------------------------
# Alphabet validation
# ---------------------------------------------------------------------------

def _validate_alphabet(seq_upper: str, alpha: str, name: str) -> None:
    """Reject sequences outside the declared alphabet. Catches amino acids."""
    legal = _LEGAL_DNA if alpha == "DNA" else _LEGAL_RNA
    for ch in seq_upper:
        if ch not in legal:
            raise ValueError(
                f"Record {name!r}: illegal symbol {ch!r} for ALPHA={alpha}. "
                f"FASTB v3 is nucleotide-only; protein sequences are not supported. "
                f"Legal characters: {sorted(legal)}"
            )


def _detect_probable_protein(seq_upper: str, name: str, threshold: float = 0.05) -> None:
    """Heuristic check for amino-acid sequences that happen to use only letters
    also valid in IUPAC nucleotide codes (e.g. silk fibroin's GAGAGS motif).

    Real nucleotide sequences are overwhelmingly ACGT/ACGU (+ occasional N);
    degenerate codes RYSWKMBDHV typically appear at <1% in real biological data.
    If the sequence is long enough to be meaningful and shows >5% "degenerate"
    characters, it is almost certainly a protein being misencoded.
    """
    if len(seq_upper) < 20:
        # Too short to apply the heuristic reliably. Common case: adapters,
        # primers, short test fixtures.
        return

    # Count letters that are legal in 4-bit mode but are NOT the four canonical
    # bases or N. High frequency of these is the protein signature.
    degenerate_letters = set("WSMKRYBDHV")
    deg_count = sum(1 for ch in seq_upper if ch in degenerate_letters)
    fraction = deg_count / len(seq_upper)

    if fraction > threshold:
        raise ValueError(
            f"Record {name!r}: {fraction:.1%} of residues are IUPAC-degenerate "
            f"codes (>{threshold:.0%} threshold). This is almost certainly a "
            f"protein sequence being misinterpreted as nucleotides. If this is "
            f"genuinely a heavily-degenerate nucleotide sequence, re-encode with "
            f"the --force-nucleotide flag."
        )


def _detect_mixed_tu(seq_upper: str, alpha: str, name: str) -> None:
    """Reject records with both T and U — the format encodes only one per record."""
    if "T" in seq_upper and "U" in seq_upper:
        raise ValueError(
            f"Record {name!r}: contains both T and U. FASTB v3 requires a single "
            f"nucleotide alphabet per record (ALPHA=DNA uses T, ALPHA=RNA uses U)."
        )
    if alpha == "DNA" and "U" in seq_upper:
        raise ValueError(
            f"Record {name!r}: declared ALPHA=DNA but sequence contains U. "
            f"Declare ALPHA=RNA or convert U to T before encoding."
        )
    if alpha == "RNA" and "T" in seq_upper:
        raise ValueError(
            f"Record {name!r}: declared ALPHA=RNA but sequence contains T. "
            f"Declare ALPHA=DNA or convert T to U before encoding."
        )


# ---------------------------------------------------------------------------
# Bit-packing
# ---------------------------------------------------------------------------

def _pick_encoding(seq_upper: str, alpha: str) -> int:
    """Choose 2 (if sequence is pure ACGT or ACGU) or 4 (any degenerate/gap)."""
    basic = "ACGT" if alpha == "DNA" else "ACGU"
    return 2 if all(ch in basic for ch in seq_upper) else 4


def _pack2(seq_upper: str) -> bytes:
    """Pack to 2 bits/base, MSB-first, last byte zero-padded."""
    n = len(seq_upper)
    out = bytearray((n + 3) // 4)
    for i, ch in enumerate(seq_upper):
        code = _PACK2[ch]
        bi = i // 4
        shift = (3 - (i % 4)) * 2
        out[bi] |= code << shift
    return bytes(out)


def _unpack2(payload: bytes, length: int, alpha: str) -> str:
    table = _UNPACK2_DNA if alpha == "DNA" else _UNPACK2_RNA
    out = []
    for i in range(length):
        bi = i // 4
        shift = (3 - (i % 4)) * 2
        code = (payload[bi] >> shift) & 0b11
        out.append(table[code])
    return "".join(out)


def _pack4(seq_upper: str) -> bytes:
    """Pack to 4 bits/base, high nibble first, last byte zero-padded."""
    n = len(seq_upper)
    out = bytearray((n + 1) // 2)
    for i, ch in enumerate(seq_upper):
        nib = _PACK4[ch]
        bi = i // 2
        if i % 2 == 0:
            out[bi] = nib << 4
        else:
            out[bi] |= nib
    return bytes(out)


def _unpack4(payload: bytes, length: int, alpha: str) -> str:
    table = _UNPACK4_DNA if alpha == "DNA" else _UNPACK4_RNA
    out = []
    for i in range(length):
        bi = i // 2
        nib = (payload[bi] >> 4) if (i % 2 == 0) else (payload[bi] & 0xF)
        out.append(table[nib])
    return "".join(out)


# ---------------------------------------------------------------------------
# Soft-masking (lowercase intervals)
# ---------------------------------------------------------------------------

def _extract_mask_runs(seq: str) -> List[Tuple[int, int]]:
    """Find (start, length) of contiguous lowercase runs."""
    runs = []
    i = 0
    n = len(seq)
    while i < n:
        if seq[i].islower():
            j = i + 1
            while j < n and seq[j].islower():
                j += 1
            runs.append((i, j - i))
            i = j
        else:
            i += 1
    return runs


def _apply_mask_runs(seq: str, runs: List[Tuple[int, int]]) -> str:
    if not runs:
        return seq
    chars = list(seq)
    for start, length in runs:
        for k in range(start, min(len(chars), start + length)):
            chars[k] = chars[k].lower()
    return "".join(chars)


def _encode_mask(runs: List[Tuple[int, int]]) -> str:
    """Hex-encoded (start:length) pairs, comma-separated. Variable-width; no 4 Gb ceiling."""
    if not runs:
        return ""
    return ",".join(f"{s:x}:{l:x}" for s, l in runs)


def _decode_mask(s: str) -> List[Tuple[int, int]]:
    if not s:
        return []
    return [(int(a, 16), int(b, 16)) for a, b in (p.split(":") for p in s.split(","))]


# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------

def write_file(records: List[Record], out, file_comment: Optional[str] = None) -> None:
    """Write records to `out` (a binary file-like).

    Args:
        records: iterable of Record objects.
        out: writable binary file-like supporting .write() and .tell().
        file_comment: optional free-text comment block, written as '#' lines
            after the magic line.
    """
    out.write(MAGIC_LINE)
    if file_comment:
        for line in file_comment.splitlines():
            out.write(b"# ")
            out.write(line.encode("utf-8"))
            out.write(b"\n")

    index_entries: List[Tuple[str, int, int]] = []  # (name, byte_offset, length)

    for rec in records:
        seq_upper = rec.sequence.upper()

        _validate_alphabet(seq_upper, rec.alpha, rec.name)
        _detect_mixed_tu(seq_upper, rec.alpha, rec.name)
        _detect_probable_protein(seq_upper, rec.name)

        enc = _pick_encoding(seq_upper, rec.alpha)
        payload = _pack2(seq_upper) if enc == 2 else _pack4(seq_upper)

        runs = rec.mask if rec.mask else _extract_mask_runs(rec.sequence)
        crc = zlib.crc32(payload) & 0xFFFFFFFF

        if rec.comment:
            comment_safe = rec.comment.replace("\n", " ").replace("\t", " ")
            out.write(b"# ")
            out.write(comment_safe.encode("utf-8"))
            out.write(b"\n")

        header_offset = out.tell()
        index_entries.append((rec.name, header_offset, len(seq_upper)))

        mask_field = _encode_mask(runs)
        header = (
            f">{rec.name}\t"
            f"LEN={len(seq_upper)}\t"
            f"ENC={enc}\t"
            f"ALPHA={rec.alpha}\t"
            f"CRC={crc:08x}\t"
            f"BYTES={len(payload)}"
        )
        if mask_field:
            header += f"\tMASK={mask_field}"
        header += "\n"
        out.write(header.encode("utf-8"))
        out.write(payload)
        out.write(b"\n")

    index_start = out.tell()
    for name, off, length in index_entries:
        out.write(f"{name}\t{off}\t{length}\n".encode("utf-8"))
    out.write(
        f"##INDEX\tENTRIES={len(index_entries)}\tOFFSETS_START={index_start}\n".encode("utf-8")
    )
    # Streaming writers leave FILE_CRC as 00000000; post-write tool can patch it in.
    out.write(b"##END\tFILE_CRC=00000000\n")


# ---------------------------------------------------------------------------
# Read
# ---------------------------------------------------------------------------

def _parse_header_line(line: bytes) -> Dict[str, str]:
    if not line.startswith(b">") or not line.endswith(b"\n"):
        raise ValueError(f"Malformed record header: {line!r}")
    body = line[1:-1].decode("utf-8")
    parts = body.split("\t")
    fields = {"NAME": parts[0]}
    for p in parts[1:]:
        if "=" not in p:
            raise ValueError(f"Malformed header field (no '='): {p!r}")
        k, v = p.split("=", 1)
        fields[k] = v
    return fields


def read_file(src) -> Iterator[Record]:
    """Iterate records from `src` (binary file-like supporting readline/read)."""
    magic = src.readline()
    if magic != MAGIC_LINE:
        raise ValueError(
            f"Not a FASTB v3 file. Expected magic {MAGIC_LINE!r}, got {magic!r}"
        )

    pending_comment: Optional[str] = None

    while True:
        line = src.readline()
        if not line:
            return  # clean EOF

        if line.startswith(b"#"):
            if line.startswith(b"##INDEX") or line.startswith(b"##END"):
                return
            # Record-scoped comment; attach to next record.
            pending_comment = line[1:].strip().decode("utf-8", errors="replace")
            continue

        if not line.startswith(b">"):
            if line.strip() == b"":
                continue
            # Index entry lines (NAME\toffset\tlength) — first one means we've
            # passed the last record. Stop.
            if b"\t" in line:
                return
            raise ValueError(f"Unexpected line outside record: {line!r}")

        fields = _parse_header_line(line)
        length = int(fields["LEN"])
        enc = int(fields["ENC"])
        alpha = fields.get("ALPHA")
        if alpha is None:
            raise ValueError(
                f"Record {fields['NAME']!r}: missing ALPHA field in header."
            )
        if alpha not in ("DNA", "RNA"):
            raise ValueError(
                f"Record {fields['NAME']!r}: ALPHA={alpha} is not supported by "
                f"this fastb version. Known values: DNA, RNA. "
                f"(AA is reserved for a future release.)"
            )
        expected_crc = int(fields["CRC"], 16)
        byte_len = int(fields["BYTES"])
        mask = _decode_mask(fields.get("MASK", ""))

        payload = src.read(byte_len)
        if len(payload) != byte_len:
            raise ValueError(
                f"Truncated payload for {fields['NAME']!r}: "
                f"expected {byte_len} bytes, got {len(payload)}"
            )

        actual_crc = zlib.crc32(payload) & 0xFFFFFFFF
        if actual_crc != expected_crc:
            raise ValueError(
                f"CRC mismatch for {fields['NAME']!r}: "
                f"header says {expected_crc:08x}, payload computes {actual_crc:08x}"
            )

        term = src.read(1)
        if term != b"\n":
            raise ValueError(
                f"Missing record terminator after {fields['NAME']!r} "
                f"(got {term!r}, expected b'\\n')"
            )

        if enc == 2:
            seq = _unpack2(payload, length, alpha)
        elif enc == 4:
            seq = _unpack4(payload, length, alpha)
        else:
            raise ValueError(f"Unknown ENC={enc} for {fields['NAME']!r}")

        if mask:
            seq = _apply_mask_runs(seq, mask)

        yield Record(
            name=fields["NAME"],
            sequence=seq,
            alpha=alpha,
            mask=mask,
            comment=pending_comment,
        )
        pending_comment = None


# ---------------------------------------------------------------------------
# Random access via footer index
# ---------------------------------------------------------------------------

def read_index(path: str) -> Dict[str, Tuple[int, int]]:
    """Read the footer index. Returns {name: (header_offset, length_in_bases)}."""
    with open(path, "rb") as f:
        # Seek to end, walk backwards looking for ##INDEX line.
        f.seek(0, 2)
        file_size = f.tell()
        # Read the last 4 KB — enough for any reasonable index trailer.
        back = min(file_size, 4096)
        f.seek(file_size - back)
        tail = f.read(back)

        # Find the ##INDEX line
        idx_marker = tail.rfind(b"\n##INDEX\t")
        if idx_marker < 0:
            raise ValueError("No ##INDEX trailer found")
        idx_line_start = file_size - back + idx_marker + 1
        f.seek(idx_line_start)
        idx_line = f.readline().decode("utf-8")

        # Parse OFFSETS_START from "##INDEX\tENTRIES=n\tOFFSETS_START=pos\n"
        parts = dict(
            p.split("=", 1) for p in idx_line.strip().split("\t")[1:]
        )
        offsets_start = int(parts["OFFSETS_START"])
        entries = int(parts["ENTRIES"])

        # Read the index entries
        f.seek(offsets_start)
        result = {}
        for _ in range(entries):
            entry = f.readline().decode("utf-8").strip()
            name, off, length = entry.split("\t")
            result[name] = (int(off), int(length))
        return result


def read_record(path: str, name: str) -> Record:
    """Read a single record by name, using the footer index for O(1) seek."""
    idx = read_index(path)
    if name not in idx:
        raise KeyError(f"Record {name!r} not in file")
    header_offset, _ = idx[name]
    with open(path, "rb") as f:
        f.seek(header_offset)
        line = f.readline()
        fields = _parse_header_line(line)
        length = int(fields["LEN"])
        enc = int(fields["ENC"])
        alpha = fields["ALPHA"]
        expected_crc = int(fields["CRC"], 16)
        byte_len = int(fields["BYTES"])
        mask = _decode_mask(fields.get("MASK", ""))

        payload = f.read(byte_len)
        if zlib.crc32(payload) & 0xFFFFFFFF != expected_crc:
            raise ValueError(f"CRC mismatch for {name!r}")

        if enc == 2:
            seq = _unpack2(payload, length, alpha)
        else:
            seq = _unpack4(payload, length, alpha)
        if mask:
            seq = _apply_mask_runs(seq, mask)

        return Record(name=name, sequence=seq, alpha=alpha, mask=mask)


# ---------------------------------------------------------------------------
# Demo / smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    recs = [
        Record("chr1", "ACGTACGTACGT", alpha="DNA", comment="Simple test record"),
        Record("chr2", "ACGTacgtNNNN", alpha="DNA", comment="Has lowercase and N"),
        Record("mito", "ACGUACGU", alpha="RNA"),
    ]

    buf = io.BytesIO()
    write_file(recs, buf, file_comment="FASTB v3 demo file\nGenerated by reference impl")
    data = buf.getvalue()

    print("=== File size ===")
    total_bases = sum(len(r.sequence) for r in recs)
    print(f"{len(data)} bytes for {total_bases} bases")

    print("\n=== Full file (as a text editor would display) ===")
    print(data.decode("utf-8", errors="replace"))

    print("=== Round-trip test ===")
    buf.seek(0)
    decoded = list(read_file(buf))
    for orig, got in zip(recs, decoded):
        ok = (
            orig.sequence == got.sequence
            and orig.alpha == got.alpha
            and orig.name == got.name
        )
        print(f"  {got.name}: {'OK' if ok else 'FAIL'}  seq={got.sequence!r}")

    # Persist to disk and verify random access
    import tempfile, os
    out_path = os.path.join(tempfile.gettempdir(), "demo.fastb")
    with open(out_path, "wb") as f:
        f.write(data)
    print(f"\nWrote {out_path} ({len(data)} bytes)")

    print("\n=== Random access by name ===")
    idx = read_index(out_path)
    print(f"Index: {idx}")
    rec = read_record(out_path, "chr2")
    print(f"Direct seek to 'chr2': {rec.sequence!r}")
