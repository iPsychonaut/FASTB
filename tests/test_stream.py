"""Chunked decode: boundaries, wrap widths, threads, and canonical run lists."""

import io
import subprocess
import sys

import numpy as np
import pytest

import fastb
from fastb import core, io as fio
from tests.test_cli import run


def _seq(n=1000, seed=0):
    """Random ACGT with N runs and lowercase runs placed to straddle chunk edges."""
    rng = np.random.default_rng(seed)
    arr = rng.choice(np.frombuffer(b"ACGT", np.uint8), n)
    for s in range(37, n - 20, 97):
        arr[s:s + 11] = ord("N")
    for s in range(5, n - 30, 61):
        arr[s:s + 23] |= 0x20
    return arr.tobytes().decode()


@pytest.fixture()
def small_chunks(monkeypatch):
    monkeypatch.setattr(core, "CHUNK_BASES", 64)
    monkeypatch.setattr(fio, "CHUNK_BASES", 64)


def _expected(records, width):
    out = b""
    for name, seq in records:
        out += b">" + name.encode() + b"\n" + fio.wrap_bytes(seq, width)
    return out


@pytest.mark.parametrize("width", [0, 7, 60, 80])
def test_stream_fasta_matches_wrap(tmp_path, small_chunks, width):
    records = [("a", _seq(1000)), ("b", "ACGT"), ("c", _seq(333, 1) + "W-")]
    path = str(tmp_path / "s.fastb")
    fastb.write_records(path, records)
    buf = io.BytesIO()
    assert fio.stream_fasta(path, buf, width) == 3
    assert buf.getvalue() == _expected(records, width)


def test_chunked_encoder_bytes_are_canonical(tmp_path, monkeypatch):
    """Same bytes whatever CHUNK_BASES is, so runs crossing chunk edges are merged."""
    seq = _seq(1000)
    a = fastb.encode_sequence(seq, "D")
    monkeypatch.setattr(core, "CHUNK_BASES", 64)
    b = fastb.encode_sequence(seq, "D")
    assert a == b
    assert fastb.decode_sequence(b[1], len(seq), b[0], "D",
                                 core._decode_runs(b[2]), core._decode_runs(b[3])) == seq


def test_threads_match_serial(tmp_path):
    records = [(f"r{i}", _seq(5000, i)) for i in range(6)]
    path = str(tmp_path / "t.fastb")
    fastb.write_records(path, records)
    serial = run("cat", path, "-w", "60").stdout
    par = run("cat", path, "-w", "60", "-p", "3")
    assert par.returncode == 0, par.stderr
    assert par.stdout == serial == _expected(records, 60).decode().replace("\n", "\n")


def test_stream_detects_crc_error(tmp_path):
    path = str(tmp_path / "c.fastb")
    fastb.write_records(path, [("a", "ACGT" * 50)])
    data = bytearray(open(path, "rb").read())
    hdr_end = data.index(b"\n", data.index(b">")) + 1
    data[hdr_end] ^= 0xFF
    open(path, "wb").write(data)
    with pytest.raises(ValueError, match="CRC"):
        fio.stream_fasta(path, io.BytesIO(), 80)
