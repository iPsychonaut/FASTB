"""Byte-exact conformance: each fixture FASTA must encode to its .fastb twice over."""

import glob
import io
import os

import pytest

import fastb
from fastb.io import iter_fasta

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")
NAMES = sorted(os.path.basename(p)[:-6] for p in glob.glob(os.path.join(FIXTURES, "*.fasta")))


def _records(name):
    recs = []
    for rec_name, seq in iter_fasta(os.path.join(FIXTURES, name + ".fasta")):
        alpha = "R" if "U" in seq.upper() else "D"
        recs.append(fastb.Record(rec_name, seq, alpha=alpha))
    return recs


@pytest.mark.parametrize("name", NAMES)
def test_encode_is_byte_exact(name):
    buf = io.BytesIO()
    fastb.write_file(_records(name), buf)
    with open(os.path.join(FIXTURES, name + ".fastb"), "rb") as f:
        assert buf.getvalue() == f.read()


@pytest.mark.parametrize("name", NAMES)
def test_decode_matches_fasta(name):
    expected = [(r.name, r.sequence) for r in _records(name)]
    with open(os.path.join(FIXTURES, name + ".fastb"), "rb") as f:
        got = [(r.name, r.sequence) for r in fastb.read_file(f)]
    assert got == expected


def test_n_runs_keep_2bit():
    """A record with N gaps must stay at ENC=2 and restore N (and n) on decode."""
    seq = "ACGT" * 10 + "N" * 20 + "acgtnnnn" + "ACGT" * 10
    enc, payload, nruns, mask, iupac = fastb.encode_sequence(seq, "D")
    assert enc == 2 and len(payload) == -(-len(seq) // 4)
    assert nruns and mask
    assert fastb.decode_sequence(payload, len(seq), 2, "D",
                                 fastb.core._decode_runs(nruns),
                                 fastb.core._decode_runs(mask)) == seq


def test_unknown_header_key_is_ignored(tmp_path):
    path = tmp_path / "x.fastb"
    buf = io.BytesIO()
    fastb.write_file([fastb.Record("s", "ACGTACGT")], buf)
    data = buf.getvalue().replace(b"\tBYTES=2", b"\tBYTES=2\tFUTURE=1", 1)
    [rec] = fastb.read_file(io.BytesIO(data))
    assert rec.sequence == "ACGTACGT"


def test_read_record_by_position(tmp_path):
    path = str(tmp_path / "k.fastb")
    with open(path, "wb") as f:
        fastb.write_file([fastb.Record("a", "AAAA"), fastb.Record("b", "CCCC"),
                          fastb.Record("c", "GGGG")], f)
    assert fastb.read_record(path, 1).sequence == "CCCC"
    assert fastb.read_record(path, "c").sequence == "GGGG"
    with pytest.raises(KeyError):
        fastb.read_record(path, 3)
