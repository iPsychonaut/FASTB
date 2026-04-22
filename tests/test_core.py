"""Tests for fastb.core: round-trips, random access, and error paths."""

import io
import zlib

import pytest

from fastb.core import (
    MAGIC_LINE,
    Record,
    read_file,
    read_index,
    read_record,
    write_file,
)


def _make_file(*records, comment=None) -> bytes:
    buf = io.BytesIO()
    write_file(list(records), buf, file_comment=comment)
    return buf.getvalue()


def _roundtrip(*records) -> list[Record]:
    data = _make_file(*records)
    return list(read_file(io.BytesIO(data)))


# ---------------------------------------------------------------------------
# Round-trips
# ---------------------------------------------------------------------------

def test_roundtrip_simple_dna():
    rec = Record("seq1", "ACGTACGTACGT", alpha="D")
    [got] = _roundtrip(rec)
    assert got.name == "seq1"
    assert got.sequence == "ACGTACGTACGT"
    assert got.alpha == "D"


def test_roundtrip_mixed_case_mask():
    """Lowercase regions must survive the round-trip as lowercase."""
    rec = Record("seq2", "ACGTacgtNNNN", alpha="D")
    [got] = _roundtrip(rec)
    assert got.sequence == "ACGTacgtNNNN"
    assert got.alpha == "D"


def test_roundtrip_rna():
    rec = Record("mito", "ACGUACGUACGU", alpha="R")
    [got] = _roundtrip(rec)
    assert got.sequence == "ACGUACGUACGU"
    assert got.alpha == "R"


def test_roundtrip_degenerates_enc4():
    """Sequences with IUPAC degenerates or gaps use ENC=4."""
    seq = "ACGTWSMKRYN-"
    rec = Record("degen", seq, alpha="D")
    [got] = _roundtrip(rec)
    assert got.sequence == seq
    assert got.alpha == "D"


# ---------------------------------------------------------------------------
# Random access
# ---------------------------------------------------------------------------

def test_random_access_middle_record(tmp_path):
    records = [
        Record("first", "AAAA", alpha="D"),
        Record("middle", "CCCCCCCC", alpha="D"),
        Record("last", "GGGGGGGG", alpha="D"),
    ]
    path = str(tmp_path / "ra.fastb")
    with open(path, "wb") as f:
        write_file(records, f)

    got = read_record(path, "middle")
    assert got.name == "middle"
    assert got.sequence == "CCCCCCCC"


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------

def test_crc_corruption_raises():
    """Flipping a payload byte must cause a CRC-mentioning ValueError."""
    rec = Record("seq", "ACGTACGT", alpha="D")
    data = bytearray(_make_file(rec))

    # Find payload start: first '>' header, skip to end of that line, flip first byte
    idx = data.index(ord(">"))
    hdr_end = data.index(ord("\n"), idx) + 1
    data[hdr_end] ^= 0xFF

    with pytest.raises(ValueError, match="CRC"):
        list(read_file(io.BytesIO(bytes(data))))


def test_truncation_raises():
    """Cutting the payload in half must cause a truncation-mentioning ValueError."""
    rec = Record("seq", "ACGTACGTACGTACGTACGT", alpha="D")
    data = _make_file(rec)

    # Keep magic + header, cut payload in half
    idx = data.index(ord(">"))
    hdr_end = data.index(ord("\n"), idx) + 1
    cut = hdr_end + 1  # keep only first payload byte of what should be 5 bytes

    with pytest.raises(ValueError, match="[Tt]runcat"):
        list(read_file(io.BytesIO(data[:cut])))


def test_unknown_alpha_value_raises():
    """A record with ALPHA=P must raise a clear ValueError mentioning the reserved value."""
    payload = bytes([0b00000000])  # 4 A's in 2-bit
    crc = zlib.crc32(payload) & 0xFFFFFFFF
    hdr = f">seq\tLEN=4\tENC=2\tALPHA=P\tCRC={crc:08x}\tBYTES=1\n"
    data = MAGIC_LINE + hdr.encode() + payload + b"\n"

    with pytest.raises(ValueError, match="P"):
        list(read_file(io.BytesIO(data)))
