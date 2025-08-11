#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTB v1 (tetrabin sentinel) and v2 (TLV container) encoder/decoder.

This module provides:
  - v1:
      * `fastb_transform`: encode/decode legacy FASTB records (sentinel-based)
      * `save_bitarray_to_file`: save aggregated FASTB v1 bitarray as binary
      * `read_fastb_file`: read a FASTB v1 file and iterate decoded records
  - v2 (preferred):
      * `fastb2_encode_records`: encode records into TLV container with CRC32
      * `fastb2_decode_stream`: decode TLV container bytes into records
      * `save_bytes_to_file`: save v2 bytes to disk
      * `read_fastb_auto`: read a file and auto-detect v2 (preferred) or v1

Record structure (v1):
  [description_bytes][type_spacer][sequence_bits][record_terminator]
Markers (v1):
  type_spacer:
    DNA -> "000011110000"
    RNA -> "111100001111"
  record_terminator:
    next_spacer -> "000000001111"

Container structure (v2, little-endian):
  File header:
    0-3   Magic:  b"FAST"
    4     Version: 0x02
    5     Flags:   bit0=has_file_index (reserved)
    6-7   Reserved: 0
    8-11  RecordCount (uint32)

  Then for each record:
    TLV blocks (Type:uint8, Len:uint32, Payload)
      0x01 DESC_UTF8   : Len=bytes, UTF-8 description
      0x02 NUC_TYPE    : Len=1, payload 'D' or 'R'
      0x03 SEQ_4B      : Len=nibbles, payload packed (two bases/byte, high nibble first)
      0x04 NOTE_UTF8   : optional notes (UTF-8)
      0x05 SOFTMASK_RLE: optional softmask runs, payload is repeated <uint32 start, uint32 length>
      0xFE END_RECORD  : optional (Len=0)
      0xFF LEGACY_FASTB: optional v1 payload (bytes) for archival
    CRC32 (uint32 LE) over the concatenated TLV bytes (not including CRC field)

Integrity:
  - v1: guards against in-band triplicate spacers when encoding
  - v2: per-record CRC32 for structured, collision-free boundaries

KISS:
  - v2 uses minimal TLV types, no scanning, and explicit lengths.
  - v1 code kept intact for compatibility/parity with existing datasets.

ISO notes:
  - UTC timestamps for logs
  - Clear input/output typing and error handling
  - Deterministic, documented transforms
"""

from __future__ import annotations

from bitarray import bitarray
from Bio import SeqIO
import pandas as pd
from rich.logging import RichHandler
import hashlib
import logging
import os
import time
import io
import struct
import zlib
import argparse
from typing import List, Tuple, Optional


# -----------------------------------------------------------------------------
# Logging (ISO-8601 UTC timestamps, color output via Rich)
# -----------------------------------------------------------------------------
_LOG_LEVEL = os.getenv("FASTB_LOGLEVEL", "INFO").upper()


class _UTCFormatter(logging.Formatter):
    """Logging formatter with ISO-8601 UTC timestamps."""

    converter = time.gmtime  # UTC

    def formatTime(self, record, datefmt=None):
        return time.strftime("%Y-%m-%dT%H:%M:%SZ", self.converter(record.created))


_formatter = _UTCFormatter("%(asctime)s | %(levelname)-8s | %(name)s | %(message)s")

logging.basicConfig(
    level=getattr(logging, _LOG_LEVEL, logging.INFO),
    handlers=[RichHandler(rich_tracebacks=True, markup=False, show_level=False, show_time=False, show_path=False)],
    force=True,
)
for _h in logging.getLogger().handlers:
    _h.setFormatter(_formatter)

logger = logging.getLogger("FASTB")


# =============================================================================
# v1 (legacy) tetrabin encode/decode
# =============================================================================
def fastb_transform(input_sequence, mode, input_description=None):
    """
    Encode or decode a nucleotide record using FASTB v1 tetrabin encoding.

    Args:
        input_sequence (str | bitarray):
            If mode == "encode": a nucleotide string (DNA or RNA).
            If mode == "decode": a bitarray representing one encoded record
            ending with the record terminator.

        mode (str):
            Either "encode" or "decode".

        input_description (str | None):
            Used only when mode == "encode": header/description text to prepend
            as raw UTF-8 bytes before the type spacer.

    Returns:
        tuple | bitarray:
            If mode == "encode":
                bitarray
                    Aggregated record:
                    [description_bytes][type_spacer][sequence_bits][record_terminator]
            If mode == "decode":
                (decoded_sequence: str, nucleotide_type: str, decoded_description: str)

    Raises:
        ValueError:
            - If mode is not "encode" or "decode".
            - If an invalid nucleotide character is encountered.
            - If encoded sequence contains a triplicate FASTB marker substring
              (three consecutive copies of DNA/RNA spacer or NEXT terminator).
            - If no valid spacer is found (decode).
            - If an invalid tetrabin block is encountered (decode).
    """
    tetrabin_scheme = {
        "U": "0100", "T": "0100", "A": "1000", "C": "0010", "G": "0001",
        "R": "1001", "Y": "0110", "K": "0101", "M": "1010", "S": "0011",
        "W": "1100", "B": "0111", "D": "1101", "H": "1110", "V": "1011",
        "N": "1111", " ": "0000"
    }

    spacer = {"DNA": "000011110000", "RNA": "111100001111"}
    next_spacer = "000000001111"
    bit_next = bitarray(next_spacer)

    inverse_scheme = {v: k for k, v in tetrabin_scheme.items()}
    spacer_reverse = {v: k for k, v in spacer.items()}

    if mode == "encode":
        cleaned_sequence = input_sequence.upper().replace("\n", "")
        logger.debug("Encoding (v1) sequence length=%d", len(cleaned_sequence))

        nucleotide_type = "RNA" if ("u" in cleaned_sequence or "U" in cleaned_sequence) else "DNA"
        bit_spacer = bitarray(spacer.get(nucleotide_type))

        try:
            encoded_sequence = "".join(tetrabin_scheme[n] for n in cleaned_sequence)
        except KeyError:
            logger.error("Invalid nucleotide detected during v1 encoding")
            raise ValueError(f"Invalid nucleotide character in sequence: {cleaned_sequence}")

        # Guard: triplicate marker detection
        triplicate_markers = {
            "DNA": spacer["DNA"] * 3,
            "RNA": spacer["RNA"] * 3,
            "NEXT": next_spacer * 3,
        }
        for label, pattern in triplicate_markers.items():
            pos = encoded_sequence.find(pattern)
            if pos != -1:
                nt_index = pos // 4
                logger.error(
                    "(v1) Triplicate FASTB marker '%s' at bit %d (nt %d); aborting.",
                    label, pos, nt_index
                )
                raise ValueError(
                    "Encoded sequence contains a triplicate FASTB marker '%s' at "
                    "bit offset %d (nt index %d); aborting." % (label, pos, nt_index)
                )

        bit_sequence = bitarray(encoded_sequence)

        bit_description = bitarray()
        bit_description.frombytes((input_description or "").encode("utf-8"))
        bit_complete = bit_description
        bit_complete.extend(bit_spacer)
        bit_complete.extend(bit_sequence)
        bit_complete.extend(bit_next)
        return bit_complete

    elif mode == "decode":
        str_sequence = input_sequence.to01() if isinstance(input_sequence, bitarray) else input_sequence
        logger.debug("Decoding (v1) record bits=%d", len(str_sequence))

        for spc in spacer_reverse:
            if spc in str_sequence:
                nucleotide_type = spacer_reverse[spc]
                prefix, coded = str_sequence.split(spc, 1)
                coded, _line_break = coded.split(next_spacer, 1)
                break
        else:
            logger.error("(v1) No valid type spacer found during decode")
            raise ValueError(f"No valid spacer found in encoded input: {str_sequence}.")

        decoded_description = (
            "".join(chr(int(prefix[i:i + 8], 2)) for i in range(0, len(prefix), 8))
            if prefix else ""
        )

        try:
            decoded_seq = "".join(
                inverse_scheme[coded[i:i + 4]]
                for i in range(0, len(coded), 4)
            )
        except KeyError:
            logger.error("(v1) Invalid tetrabin block during decode")
            raise ValueError("Invalid tetrabin block in encoded sequence.")

        if nucleotide_type == "RNA":
            decoded_seq = decoded_seq.replace("T", "U")

        return decoded_seq, nucleotide_type, decoded_description

    else:
        logger.error("Invalid mode '%s'", mode)
        raise ValueError("Mode must be 'encode' or 'decode'.")


def save_bitarray_to_file(bit_data: bitarray, output_path: str) -> None:
    """Save a bitarray object as a binary file (v1 utility)."""
    if not isinstance(bit_data, bitarray):
        logger.error("save_bitarray_to_file received non-bitarray input")
        raise TypeError("Input must be a bitarray object.")
    try:
        with open(output_path, 'wb') as file:
            file.write(bit_data.tobytes())
        logger.info("FASTB file written: %s", output_path)
    except IOError as e:
        logger.error("Failed to write FASTB file: %s", output_path)
        raise IOError(f"Failed to write to file: {output_path}") from e


def read_fastb_file(fastb_path: str) -> List[Tuple[str, str, str]]:
    """Read a FASTB v1 file and iterate over records, decoding them to verify integrity."""
    next_spacer = "000000001111"
    bit_data = bitarray()
    try:
        with open(fastb_path, 'rb') as file:
            bit_data.fromfile(file)
        logger.info("FASTB file read: %s", fastb_path)
    except IOError as e:
        logger.error("Failed to read FASTB file: %s", fastb_path)
        raise IOError(f"Failed to read FASTB file: {fastb_path}") from e

    str_data = bit_data.to01()
    records: List[Tuple[str, str, str]] = []
    chunks = str_data.split(next_spacer)
    logger.debug("(v1) Split FASTB stream into %d chunks", len(chunks))

    for i, chunk in enumerate(chunks, start=1):
        if not chunk.strip("0"):
            logger.debug("(v1) Chunk %d: skipped (all zeros / empty)", i)
            continue

        chunk_bits = bitarray(chunk + next_spacer)
        try:
            seq_str, nuc_type, desc = fastb_transform(chunk_bits, "decode")
            records.append((desc, seq_str, nuc_type))
        except ValueError as e:
            logger.error("(v1) Decode error in chunk %d: %s", i, e)

    logger.info("(v1) Decoded %d FASTB record(s) from %s", len(records), fastb_path)
    return records


# =============================================================================
# v2 (preferred) TLV container encode/decode
# =============================================================================

# File header
_FASTB2_MAGIC = b"FAST"
_FASTB2_VER = 0x02

# TLV Types
_DESC_UTF8 = 0x01
_NUC_TYPE = 0x02   # payload: b"D" or b"R"
_SEQ_4B = 0x03     # length is number of nibbles, payload packed (two bases/byte)
_NOTE_UTF8 = 0x04
_SOFTMASK_RLE = 0x05  # payload: repeated <uint32 start, uint32 length> little-endian
_END_RECORD = 0xFE
_LEGACY_FASTB = 0xFF  # optional archival v1 payload

# Tetrabin mapping
_TETRA_NIB = {
    "U": 0x4, "T": 0x4, "A": 0x8, "C": 0x2, "G": 0x1,
    "R": 0x9, "Y": 0x6, "K": 0x5, "M": 0xA, "S": 0x3,
    "W": 0xC, "B": 0x7, "D": 0xD, "H": 0xE, "V": 0xB,
    "N": 0xF, " ": 0x0,
}
_INV_NIB = {v: k for k, v in _TETRA_NIB.items()}


def _find_softmask_runs(seq: str) -> List[Tuple[int, int]]:
    """Return [(start, length), ...] for contiguous lowercase runs."""
    runs: List[Tuple[int, int]] = []
    i, n = 0, len(seq)
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


def _apply_softmask_runs(seq: str, runs: List[Tuple[int, int]]) -> str:
    """Lowercase seq positions per runs; returns new string."""
    if not runs:
        return seq
    s = list(seq)
    L = len(s)
    for start, length in runs:
        end = min(L, start + length)
        for k in range(start, end):
            s[k] = s[k].lower()
    return "".join(s)


def _pack_seq_4b(seq: str) -> Tuple[bytes, int]:
    """
    Pack a nucleotide sequence to 4-bit nibbles (two bases per byte).

    Args:
        seq (str): Nucleotide string.

    Returns:
        (packed_bytes, nibble_count)
    """
    s_raw = (seq or "").replace("\n", "").replace(" ", "")
    s = s_raw.upper()  # mapping only; case preserved separately via SOFTMASK_RLE
    try:
        nibbles = [_TETRA_NIB[ch] for ch in s]
    except KeyError as e:
        raise ValueError(f"Invalid nucleotide '{e.args[0]}' in sequence.") from None

    nib_count = len(nibbles)
    out = bytearray((nib_count + 1) // 2)
    for i, nib in enumerate(nibbles):
        bi = i // 2
        if i % 2 == 0:
            out[bi] = (nib << 4)
        else:
            out[bi] |= nib
    return bytes(out), nib_count


def _unpack_seq_4b(buf: bytes, nib_count: int, nuc_type: str) -> str:
    """
    Unpack 4-bit packed bases into a sequence string.

    Args:
        buf (bytes): Packed bytes.
        nib_count (int): Number of nibbles to read.
        nuc_type (str): 'D' for DNA, 'R' for RNA.

    Returns:
        str: Decoded sequence in canonical letters (RNA uses 'U').
    """
    bases = []
    for i in range(nib_count):
        b = buf[i // 2]
        nib = (b >> 4) if (i % 2 == 0) else (b & 0xF)
        base = _INV_NIB.get(nib)
        if base is None:
            raise ValueError(f"Invalid tetrabin nibble 0x{nib:X}")
        bases.append(base)
    s = "".join(bases)
    return s.replace("T", "U") if nuc_type == "R" else s


def _write_tlv(w: io.BufferedWriter, t: int, length: int, payload: Optional[bytes] = None) -> None:
    """Write one TLV to a buffer."""
    w.write(struct.pack("<BI", t, length))
    if length and payload:
        w.write(payload)


def fastb2_encode_records(records: List[Tuple[str, str, str]], include_legacy_v1: bool = False) -> bytes:
    """
    Encode records into FASTB v2 TLV container (preferred).

    Each record stores DESC, NUC_TYPE, SEQ_4B, and (if needed) SOFTMASK_RLE.

    Args:
        records: list of (description, sequence, nucleotide_type)
                 nucleotide_type is 'DNA' or 'RNA' (case-insensitive).
        include_legacy_v1: If True, embeds a legacy v1 bitstream TLV per record.

    Returns:
        bytes: Encoded file bytes.
    """
    out = io.BytesIO()
    out.write(_FASTB2_MAGIC)
    out.write(struct.pack("<B B H I", _FASTB2_VER, 0, 0, len(records)))  # ver, flags, reserved, count

    for desc, seq, nuc in records:
        nuc_tag = "R" if str(nuc).upper().startswith("R") or ("U" in (seq or "").upper()) else "D"

        rec_buf = io.BytesIO()
        desc_b = (desc or "").encode("utf-8")
        _write_tlv(rec_buf, _DESC_UTF8, len(desc_b), desc_b)
        _write_tlv(rec_buf, _NUC_TYPE, 1, nuc_tag.encode("ascii"))

        # Sequence packing (uppercased for mapping only)
        packed, nibs = _pack_seq_4b(seq or "")
        _write_tlv(rec_buf, _SEQ_4B, nibs, packed)

        # Optional softmask runs: indices align to cleaned (no spaces/newlines) original
        seq_clean = (seq or "").replace("\n", "").replace(" ", "")
        runs = _find_softmask_runs(seq_clean)
        if runs:
            payload = b"".join(struct.pack("<II", start, length) for start, length in runs)
            _write_tlv(rec_buf, _SOFTMASK_RLE, len(payload), payload)

        if include_legacy_v1:
            v1_bits = fastb_transform(seq or "", "encode", desc or "")
            _write_tlv(rec_buf, _LEGACY_FASTB, len(v1_bits.tobytes()), v1_bits.tobytes())

        rec_bytes = rec_buf.getvalue()
        crc = zlib.crc32(rec_bytes) & 0xFFFFFFFF
        out.write(rec_bytes)
        out.write(struct.pack("<I", crc))

    return out.getvalue()


def fastb2_decode_stream(data: bytes) -> List[Tuple[str, str, str]]:
    """
    Decode a FASTB v2 TLV container byte stream into records.

    Returns:
        list[(description, sequence, 'DNA'|'RNA')]
    """
    rd = io.BytesIO(data)
    if rd.read(4) != _FASTB2_MAGIC:
        raise ValueError("Not a FASTB v2 container (magic mismatch).")
    ver, _flags, _res, count = struct.unpack("<B B H I", rd.read(8))
    if ver != _FASTB2_VER:
        raise ValueError(f"Unsupported FASTB container version: {ver}")

    records: List[Tuple[str, str, str]] = []
    for idx in range(count):
        tlv_bytes = bytearray()
        desc: str = ""
        nuc_tag: Optional[str] = None
        seq_buf: bytes = b""
        nibs: int = 0
        softmask_runs: List[Tuple[int, int]] = []

        # Read TLVs until CRC (writer emits CRC immediately after TLVs)
        while True:
            hdr = rd.read(5)
            if len(hdr) < 5:
                raise ValueError(f"Unexpected EOF in TLV header for record {idx+1}")
            t, L = struct.unpack("<BI", hdr)

            # Payload length: for SEQ_4B, payload bytes = ceil(L/2)
            pay_len = ((L + 1) // 2) if t == _SEQ_4B else L
            payload = rd.read(pay_len)
            if len(payload) < pay_len:
                raise ValueError(f"Unexpected EOF in TLV payload for record {idx+1}")

            tlv_bytes += hdr + payload

            if t == _DESC_UTF8:
                desc = payload.decode("utf-8")
            elif t == _NUC_TYPE:
                nuc_tag = payload.decode("ascii")
            elif t == _SEQ_4B:
                seq_buf = payload
                nibs = L
            elif t == _SOFTMASK_RLE:
                # parse <uint32 start, uint32 length> tuples
                if L % 8 != 0:
                    raise ValueError(f"Malformed SOFTMASK_RLE length in record {idx+1}")
                softmask_runs = [
                    struct.unpack("<II", payload[k:k + 8]) for k in range(0, L, 8)
                ]
            elif t in (_NOTE_UTF8, _LEGACY_FASTB, _END_RECORD):
                pass  # optional/ignored here
            else:
                pass  # forward-compatible: ignore unknown TLVs

            # Heuristic: try reading CRC next
            pos = rd.tell()
            crc_bytes = rd.read(4)
            if len(crc_bytes) < 4:
                rd.seek(pos)
                continue
            crc_read = struct.unpack("<I", crc_bytes)[0]
            crc_calc = zlib.crc32(bytes(tlv_bytes)) & 0xFFFFFFFF
            if crc_calc != crc_read:
                rd.seek(pos)
                continue
            # TLVs complete for this record
            break

        if nuc_tag not in ("D", "R"):
            raise ValueError(f"Record {idx+1}: missing/invalid NUC_TYPE")

        # Rebuild sequence, then apply optional softmask
        seq = _unpack_seq_4b(seq_buf, nibs, nuc_tag)
        if softmask_runs:
            seq = _apply_softmask_runs(seq, softmask_runs)

        records.append((desc, seq, "RNA" if nuc_tag == "R" else "DNA"))

    return records


def save_bytes_to_file(data: bytes, output_path: str) -> None:
    """Save raw bytes to a file (v2 utility)."""
    try:
        with open(output_path, "wb") as f:
            f.write(data)
        logger.info("FASTB v2 file written: %s", output_path)
    except IOError as e:
        logger.error("%s\tFailed to write FASTB v2 file: %s", e, output_path)
        raise


def read_fastb_auto(path: str) -> List[Tuple[str, str, str]]:
    """
    Auto-detect and decode FASTB file:
      - Prefer v2 container (magic b'FAST', ver=0x02).
      - Fallback to v1 sentinel format.
    """
    with open(path, "rb") as f:
        head = f.read(12)
        f.seek(0)
        data = f.read()

    if len(head) >= 5 and head[:4] == _FASTB2_MAGIC and head[4] == _FASTB2_VER:
        logger.info("Detected FASTB v2 container: %s", path)
        return fastb2_decode_stream(data)

    logger.info("Falling back to FASTB v1 sentinel decoding: %s", path)
    return read_fastb_file(path)


def _find_first_diff(a: str, b: str) -> int:
    m = min(len(a), len(b))
    for i in range(m):
        if a[i] != b[i]:
            return i
    return -1 if len(a) == len(b) else m


def _mismatch_reason(orig: str | None, dec: str) -> str:
    if orig is None:
        return "desc-missing-in-original"
    if len(orig) != len(dec):
        return "length-mismatch"
    if orig.replace('U', 'T') == dec.replace('U', 'T'):
        has_u = 'U' in orig
        has_t = 'T' in orig
        has_u_dec = 'U' in dec
        has_t_dec = 'T' in dec
        if has_u and has_t:
            return "mixed-T/U-in-original (lossy 0x4 nibble)"
        if has_u != has_u_dec or has_t != has_t_dec:
            return "global-T↔U-normalization"
        return "T/U-only-difference"
    return "base-content-mismatch"


def fasta_conversion(input_fasta: str) -> str:
    fasta_df = pd.DataFrame(
        [(rec.description, str(rec.seq)) for rec in SeqIO.parse(input_fasta, "fasta")],
        columns=["Description", "Sequence"]
    )
    output_fastb = input_fasta.replace(".fasta", ".fastb")

    # Pre-scan mixed T/U (heads-up only)
    mixed = []
    for desc, seq in zip(fasta_df["Description"], fasta_df["Sequence"]):
        su = seq.upper()
        if ("T" in su) and ("U" in su):
            mixed.append(desc)
    if mixed:
        logger.warning(
            "Detected %d record(s) with mixed T and U (cannot be losslessly reconstructed from 4-bit 0x4):",
            len(mixed)
        )
        for d in mixed[:10]:
            logger.warning("  mixed T/U: %s", (d[:200] + "…") if len(d) > 200 else d)
        if len(mixed) > 10:
            logger.warning("  …and %d more", len(mixed) - 10)

    # Build records list once
    records: List[Tuple[str, str, str]] = []
    for _, row in fasta_df.iterrows():
        desc = row["Description"]
        seq = row["Sequence"]
        nuc = "RNA" if ("u" in seq or "U" in seq) else "DNA"
        records.append((desc, seq, nuc))

    # Encode (v2 preferred)
    if os.getenv("FASTB_LEGACY", "0") == "1":
        logger.info("FASTB_LEGACY=1 -> emitting legacy v1 bitstream")
        all_binary_data = bitarray()
        for desc, seq, _nuc in records:
            v1_bits = fastb_transform(seq, "encode", desc)
            # sanity round-trip
            seq_str, _nuc_type, desc2 = fastb_transform(v1_bits, "decode")
            assert desc2 == desc and (seq_str == seq), "v1 round-trip mismatch"
            all_binary_data.extend(v1_bits)
        save_bitarray_to_file(all_binary_data, output_fastb)
    else:
        logger.info("Emitting FASTB v2 TLV container (preferred)")
        v2_bytes = fastb2_encode_records(records, include_legacy_v1=False)
        save_bytes_to_file(v2_bytes, output_fastb)

    # Read-back (auto-detect v2/v1) and list records
    decoded_records = read_fastb_auto(output_fastb)
    logger.info("%d record(s) read back", len(decoded_records))
    for idx, (desc, seq, nuc) in enumerate(decoded_records, start=1):
        logger.debug("Record %d | %s | %s | len=%d", idx, nuc, desc, len(seq))

    # Integrity comparison against original FASTA by Description + Sequence
    desc_to_seq = {row["Description"]: row["Sequence"] for _, row in fasta_df.iterrows()}
    seq_to_descs: dict[str, set[str]] = {}
    for _, row in fasta_df.iterrows():
        seq_to_descs.setdefault(row["Sequence"], set()).add(row["Description"])
    desc_to_sha = {d: hashlib.sha256(s.encode("utf-8")).hexdigest()
                   for d, s in desc_to_seq.items()}

    mismatch_count = 0
    seq_not_found_count = 0
    pair_mismatch_count = 0
    len_mismatch_count = 0
    hash_mismatch_count = 0

    details = []
    max_show = 10

    for idx, (desc, seq, _nuc) in enumerate(decoded_records, start=1):
        original_seq = desc_to_seq.get(desc)
        ok = (original_seq == seq)
        if not ok:
            mismatch_count += 1
            reason = _mismatch_reason(original_seq, seq)
            i = _find_first_diff(original_seq or "", seq)
            if i >= 0:
                start = max(0, i - 20)
                end = min(len(seq), i + 20)
                ctx_orig = (original_seq or "")[start:end]
                ctx_dec = seq[start:end]
            else:
                ctx_orig = ctx_dec = ""
            details.append({
                "index": idx,
                "description": desc,
                "reason": reason,
                "first_diff_at": i,
                "orig_T": (original_seq or "").count("T"),
                "orig_U": (original_seq or "").count("U"),
                "dec_T": seq.count("T"),
                "dec_U": seq.count("U"),
                "ctx_orig": ctx_orig,
                "ctx_dec": ctx_dec,
            })
        logger.debug("Round-trip equality [%d]: %s", idx, ok)

        if seq not in seq_to_descs:
            seq_not_found_count += 1
            logger.debug("Sequence presence [%d]: NOT FOUND", idx)
        else:
            logger.debug("Sequence presence [%d]: FOUND (%d description(s))", idx, len(seq_to_descs[seq]))
            pair_ok = desc in seq_to_descs[seq]
            if not pair_ok:
                pair_mismatch_count += 1
            logger.debug("Pairing check [%d]: %s", idx, pair_ok)

        if original_seq is not None:
            len_ok = (len(original_seq) == len(seq))
            if not len_ok:
                len_mismatch_count += 1
            logger.debug(
                "Length check [%d]: %s (orig=%s, dec=%d)",
                idx, len_ok, len(original_seq), len(seq)
            )

            orig_sha = desc_to_sha.get(desc)
            dec_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()
            hash_ok = (orig_sha == dec_sha)
            if not hash_ok:
                hash_mismatch_count += 1
            logger.debug("SHA256 check [%d]: %s", idx, hash_ok)

    logger.info("Integrity summary:")
    logger.info("  Description->Sequence mismatches: %d", mismatch_count)
    logger.info("  Sequence NOT FOUND in original FASTA: %d", seq_not_found_count)
    logger.info("  (desc, seq) pairing mismatches: %d", pair_mismatch_count)
    logger.info("  Length mismatches: %d", len_mismatch_count)
    logger.info("  SHA256 mismatches: %d", hash_mismatch_count)

    if any([mismatch_count, seq_not_found_count, pair_mismatch_count, len_mismatch_count, hash_mismatch_count]):
        logger.error("Integrity check: ISSUES DETECTED (see counts above).")
    else:
        logger.info("Integrity check: all decoded records match originals by desc & seq.")

    if details:
        logger.warning("---- Detailed mismatch report (showing up to %d) ----", min(max_show, len(details)))
        for d in details[:max_show]:
            logger.warning(
                "Record %d | reason=%s | first_diff=%s | T/U orig=%d/%d dec=%d/%d | desc=%s",
                d["index"], d["reason"], d["first_diff_at"],
                d["orig_T"], d["orig_U"], d["dec_T"], d["dec_U"],
                (d["description"][:120] + "…") if len(d["description"]) > 120 else d["description"]
            )
            if d["first_diff_at"] != -1:
                logger.warning("  ctx orig: %s", d["ctx_orig"])
                logger.warning("  ctx dec : %s", d["ctx_dec"])

    return output_fastb


# -----------------------------------------------------------------------------
# Script usage (v2 preferred; set FASTB_LEGACY=1 to force v1 output)
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    print("""
                          \033[0m,;\033[94m%%%%%%%%%,   \033[0m,;\033[94m%%%%%%%,   \033[0m,;\033[94m%%%%%%%, \033[0m,;\033[94m%%%%%%%%%%%, \033[0m,;\033[94m%%%%%%%%%%,
  \033[96m███████████\033[93m███████████  \033[0m;;\033[94m%%%%%%%%%% \033[0m;;\033[94m%%%%%%%%%%% \033[0m;;\033[94m%%%%%%%%% \033[0m;;\033[94m%%%%%%%%%%%% \033[0m;;\033[94m%%%%%%%%%%%
  \033[96m███████████\033[93m███████████  \033[0m;%\033[94m@@@\033[0m''''''' \033[0m;%\033[94m@@@\033[0m''';%\033[94m@@@ \033[0m;%\033[94m@@@\033[0m'''''   '''';%\033[94m@@@\033[0m'''' ;%\033[94m@@@\033[0m'''''\033[94m###
  \033[96m███████████\033[93m███████████  \033[0m;%\033[94m@@@        \033[0m;%\033[94m@@@   \033[0m;%\033[94m@@@ \033[0m;%\033[94m@@@            \033[0m;%\033[94m@@@     \033[0m;%\033[94m@@@    \033[0m;%\033[94m###
  \033[96m███████████\033[93m███████████  \033[0m;@\033[94m#########  \033[0m;@\033[94m########### \033[0m;@\033[94m#####,         \033[0m;@\033[94m###     \033[0m;@\033[94m###   \033[0m;@\033[94m###
  \033[96m███████████\033[93m███████████  \033[0m;@\033[94m#########  \033[0m;@\033[94m###########  \033[0m;@\033[94m#######,,     \033[0m;@\033[94m###     \033[0m;@\033[94m##########
  \033[92m███████████\033[91m███████████  \033[0m;@\033[94m###\033[0m;;;;;'  \033[0m;@\033[94m###\033[0m;;;;\033[94m####   \033[0m';;@\033[94m######;    \033[0m;@\033[94m###     \033[0m;@\033[94m###########
  \033[92m███████████\033[91m███████████  \033[0m;@\033[94m###        \033[0m;@\033[94m###    ####      \033[0m';;@\033[94m####    \033[0m;@\033[94m###     \033[0m;@\033[94m###\033[0m;;;;;;\033[94m###
  \033[92m███████████\033[91m███████████  \033[0m;@\033[94m###        \033[0m;@\033[94m###    ####       \033[0m';@\033[94m###'    \033[0m;@\033[94m###     \033[0m;@\033[94m###     \033[0m;@\033[94m###
  \033[92m███████████\033[91m███████████  \033[0m;#\033[94m888        \033[0m;#\033[94m888    8888 \033[0m;#\033[94m8888888888     \033[0m;#\033[94m888     \033[0m;#\033[94m888888888888
  \033[92m███████████\033[91m███████████  \033[0m;#\033[94m$$$        \033[0m;#\033[94m$$$    $$$$  \033[0m`#\033[94m$$$$$$$$      \033[0m;#\033[94m$$$     \033[0m;#\033[94m$$$$$$$$$$#
                              \033[94m╔════════════════════════════════════════════════════════════╗
                              \033[94m║\033[0m          Fast-Binary Nucleotide Encoder & Decoder          \033[94m║
                              \033[94m╚════════════════════════════════════════════════════════════╝
  
                         Curated & Maintained by Ian M Bollinger               
                              (\033[94mian.bollinger@entheome.org)\033[0m
  
                                        fastb.py
                                       version  2
    """)
    
    parser = argparse.ArgumentParser(
        description="Fast-Binary Nucleotide Encoding & Decoding"
    )
    
    # Establish Default Values for fallback
    default_input_fasta = "/run/media/EYE/Main/TESTING_SPACE/nextflow_test/Psilocybe_subtropicalis/Psilocybe_subtropicalis-Entheome/Psilocybe_subtropicalis-Entheome_final_curated.fasta" # "/var/home/EYE/Desktop/Python Programs/FASTB/OR140556.fasta"
    
    # Parse Input Arguments    
    parser.add_argument(
        "-i", "--input",
        type=str,
        default=default_input_fasta,
        help="Path to the input FASTA file"
    )
    args = parser.parse_args()
    input_fasta = args.input
    
    # Set FASTB_LEGACY=1 if you must emit legacy v1 bitstream instead. v1 Loses lower-case/Confidence nucleotide calls.
    output_fastb = fasta_conversion(input_fasta)