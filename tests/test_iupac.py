"""IUPAC exception list (format 3.1), the ENC choice rule, and reading 3.0 files."""

import glob
import io
import os

import numpy as np
import pytest

import fastb
from fastb import core, io as fio

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")
V30 = sorted(os.path.basename(p)[:-6] for p in glob.glob(os.path.join(FIXTURES, "v3.0", "*.fasta")))


def _roundtrip(seq, alpha="D"):
    enc, payload, nruns, mask, iupac = fastb.encode_sequence(seq, alpha)
    back = fastb.decode_sequence(payload, len(seq), enc, alpha, core._decode_runs(nruns),
                                 core._decode_runs(mask), core._decode_iupac(iupac))
    return enc, iupac, back


def _contig(n=2000, seed=0):
    rng = np.random.default_rng(seed)
    return rng.choice(np.frombuffer(b"ACGT", np.uint8), n).tobytes().decode()


def test_a_few_codes_stay_2bit():
    s = list(_contig())
    s[100] = "K"; s[900] = "R"; s[901] = "R"; s[1500] = "y"; s[1700:1710] = "N" * 10
    seq = "".join(s)
    enc, iupac, back = _roundtrip(seq)
    assert enc == 2
    assert iupac == "64:1:K,384:2:R,5dc:1:Y"
    assert back == seq                       # lowercase y comes back through MASK


def test_gap_runs_and_dot():
    seq = _contig(400)[:100] + "----" + _contig(400)[100:300] + ".." + "ACGT"
    enc, iupac, back = _roundtrip(seq)
    assert enc == 2 and iupac == "64:4:-,130:2:-"
    assert back == seq.replace(".", "-")     # '.' is written as gap, as in 3.0


def test_dense_record_uses_4bit():
    seq = "ACGTRYKM" * 50                     # half the bases are ambiguity codes
    enc, iupac, back = _roundtrip(seq)
    assert enc == 4 and iupac == "" and back == seq


@pytest.mark.parametrize("n,expect", [(44, 4), (46, 4), (47, 2), (400, 2)])
def test_enc_choice_rule(n, expect):
    """ENC=2 when len('\\tIUPAC=' + text) <= ceil(n/2) - ceil(n/4).

    One K at position 10 gives the text 'a:1:K': 7 + 5 = 12 bytes. 4-bit costs
    11 extra bytes at 46 bases and 12 at 47, so 47 is the first length where
    the list is not larger.
    """
    s = list("ACGT" * (n // 4) + "ACGT"[: n % 4]); s[10] = "K"
    enc, iupac, back = _roundtrip("".join(s))
    assert (12 <= -(-n // 2) - -(-n // 4)) == (expect == 2)
    assert enc == expect and back == "".join(s)


def test_rna_with_codes():
    seq = "ACGU" * 100 + "R" + "ACGU" * 100
    enc, iupac, back = _roundtrip(seq, "R")
    assert enc == 2 and iupac == "190:1:R" and back == seq


def test_bytes_do_not_depend_on_chunk_size(monkeypatch):
    s = list(_contig(1000))
    s[60:70] = "K" * 10                       # run crosses the 64-base chunk edge
    s[127] = "R"; s[128] = "R"                # run crosses the next edge
    s[300] = "W"
    seq = "".join(s)
    whole = fastb.encode_sequence(seq, "D")
    monkeypatch.setattr(core, "CHUNK_BASES", 64)
    chunked = fastb.encode_sequence(seq, "D")
    assert whole == chunked
    assert whole[4] == "3c:a:K,7f:2:R,12c:1:W"


@pytest.mark.parametrize("width,threads", [(0, 1), (60, 1), (80, 1), (60, 3)])
def test_streaming_decode_with_codes(tmp_path, monkeypatch, width, threads):
    monkeypatch.setattr(core, "CHUNK_BASES", 64)
    monkeypatch.setattr(fio, "CHUNK_BASES", 64)
    s = list(_contig(1000)); s[60:70] = "K" * 10; s[500] = "r"; s[700:720] = "N" * 20
    records = [("a", "".join(s)), ("b", _contig(300, 1))]
    path = str(tmp_path / "x.fastb")
    fastb.write_records(path, records)
    with open(path, "rb") as f:
        assert f.readline() == b"##FASTB 3.1\n"
    buf = io.BytesIO()
    if threads == 1:
        fio.stream_fasta(path, buf, width)
    else:
        monkeypatch.undo()                    # worker processes do not see monkeypatches
        fio.stream_fasta(path, buf, width, threads)
    expected = b"".join(b">" + n.encode() + b"\n" + fio.wrap_bytes(q, width) for n, q in records)
    assert buf.getvalue() == expected


@pytest.mark.parametrize("name", V30)
def test_reads_3_0_files(name):
    """Files written before the IUPAC list existed must still decode exactly."""
    path = os.path.join(FIXTURES, "v3.0", name + ".fastb")
    with open(path, "rb") as f:
        assert f.readline() == b"##FASTB 3.0\n"
    expected = list(fastb.iter_fasta(os.path.join(FIXTURES, "v3.0", name + ".fasta")))
    assert list(fastb.iter_records(path)) == expected


def test_unknown_version_is_refused():
    buf = io.BytesIO()
    fastb.write_file([fastb.Record("s", "ACGT")], buf)
    data = buf.getvalue().replace(b"##FASTB 3.1\n", b"##FASTB 3.2\n", 1)
    with pytest.raises(ValueError, match="FASTB"):
        list(fastb.read_file(io.BytesIO(data)))


def test_illegal_symbol_still_rejected():
    with pytest.raises(ValueError, match="illegal symbol"):
        fastb.encode_sequence("ACGT" * 50 + "X" + "ACGT" * 50, "D")
