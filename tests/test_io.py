"""Streaming API, wrapping, and the cat/index/info subcommands."""

import gzip
import os

import pytest

import fastb
from fastb.io import wrap_bytes
from tests.test_cli import run


@pytest.mark.parametrize("seq,width,expected", [
    ("", 80, b"\n"),
    ("ACGT", 0, b"ACGT\n"),
    ("ACGT", 4, b"ACGT\n"),
    ("ACGTA", 4, b"ACGT\nA\n"),
    ("ACGTACGT", 3, b"ACG\nTAC\nGT\n"),
])
def test_wrap_bytes(seq, width, expected):
    assert wrap_bytes(seq, width) == expected


def test_write_and_iter_records_stream(tmp_path):
    path = str(tmp_path / "s.fastb")
    n = fastb.write_records(path, ((f"r{i}", "ACGT" * (i + 1)) for i in range(5)))
    assert n == 5
    assert list(fastb.iter_records(path)) == [(f"r{i}", "ACGT" * (i + 1)) for i in range(5)]


def test_write_records_leaves_no_partial_file(tmp_path):
    path = str(tmp_path / "bad.fastb")

    def gen():
        yield "ok", "ACGT"
        yield "protein", "MKWVTFISLLLLFSSAYSRGVFRR"

    with pytest.raises(fastb.AlphabetError):
        fastb.write_records(path, gen())
    assert not os.path.exists(path) and not os.path.exists(path + ".tmp")


def test_encode_accepts_gz(tmp_path):
    fa = tmp_path / "in.fasta.gz"
    with gzip.open(fa, "wb") as f:
        f.write(b">a desc\nACGTNNNNacgt\n>b\nGGCC\n")
    result = run("encode", str(fa))
    assert result.returncode == 0, result.stderr
    out = tmp_path / "in.fastb"
    assert list(fastb.iter_records(str(out))) == [("a", "ACGTNNNNacgt"), ("b", "GGCC")]


def test_cat_matches_input(tmp_path):
    fa = tmp_path / "in.fna"
    fa.write_bytes(b">a\nACGTNNNNacgt\n>b\nGGCC\n")
    assert run("encode", str(fa)).returncode == 0
    result = run("cat", str(tmp_path / "in.fastb"), "-w", "0")
    assert result.stdout.replace("\r\n", "\n") == ">a\nACGTNNNNacgt\n>b\nGGCC\n"


def test_index_and_info(tmp_path):
    path = str(tmp_path / "x.fastb")
    fastb.write_records(path, [("a", "ACGTNNNN"), ("b", "ACGTW" + "ACGT" * 10)])
    idx = run("index", path).stdout.splitlines()
    assert idx[0].split("\t") == ["record", "name", "offset", "length"]
    assert idx[1].split("\t")[:2] == ["0", "a"] and idx[2].split("\t")[:2] == ["1", "b"]
    info = dict(line.split("\t") for line in run("info", path).stdout.splitlines())
    assert info["records"] == "2" and info["bases"] == "53"
    assert info["records_2bit"] == "1" and info["records_4bit"] == "1"
    assert info["n_runs"] == "1"
