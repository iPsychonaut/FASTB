"""Tests for the fastb CLI via subprocess."""

import os
import subprocess
import sys
import tempfile

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTB_CMD = [sys.executable, "-m", "fastb.cli"]


def run(*args, input_text=None):
    env = {**os.environ, "PYTHONPATH": REPO_ROOT}
    return subprocess.run(
        FASTB_CMD + list(args),
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        env=env,
        input=input_text,
    )


@pytest.fixture()
def small_fasta(tmp_path):
    fa = tmp_path / "test.fasta"
    fa.write_text(">chr1\nACGTACGTACGT\n>chr2\nACGTacgtNNNN\n>chr3\nACGUACGU\n")
    return str(fa)


@pytest.fixture()
def small_fastb(tmp_path, small_fasta):
    fb = str(tmp_path / "test.fastb")
    result = run("encode", small_fasta, "-o", fb)
    assert result.returncode == 0, result.stderr
    return fb


# ---------------------------------------------------------------------------
# encode
# ---------------------------------------------------------------------------

def test_encode_produces_valid_file(small_fasta, tmp_path):
    fb = str(tmp_path / "out.fastb")
    result = run("encode", small_fasta, "-o", fb)
    assert result.returncode == 0
    with open(fb, "rb") as f:
        magic = f.read(12)
    assert magic == b"##FASTB 3.1\n"


def test_encode_protein_fasta_exits_4(tmp_path):
    fa = tmp_path / "protein.fasta"
    fa.write_text(">silk\nMKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHR\n")
    fb = str(tmp_path / "protein.fastb")
    result = run("encode", str(fa), "-o", fb)
    assert result.returncode == 4
    assert "amino acid" in result.stderr.lower()
    assert not os.path.exists(fb), "partial output file must not be created"


# ---------------------------------------------------------------------------
# view
# ---------------------------------------------------------------------------

def test_view_reproduces_sequences(small_fastb):
    result = run("view", small_fastb)
    assert result.returncode == 0
    assert ">chr1" in result.stdout
    assert "ACGTACGTACGT" in result.stdout
    assert ">chr2" in result.stdout
    assert ">chr3" in result.stdout


# ---------------------------------------------------------------------------
# head
# ---------------------------------------------------------------------------

def test_head_returns_one_record(small_fastb):
    result = run("head", "-n", "1", small_fastb)
    assert result.returncode == 0
    assert result.stdout.count(">") == 1
    assert ">chr1" in result.stdout


# ---------------------------------------------------------------------------
# extract
# ---------------------------------------------------------------------------

def test_extract_returns_named_record(small_fastb):
    result = run("extract", small_fastb, "chr2")
    assert result.returncode == 0
    assert ">chr2" in result.stdout
    assert ">chr1" not in result.stdout
    assert ">chr3" not in result.stdout


def test_extract_missing_name_exits_2(small_fastb):
    result = run("extract", small_fastb, "nonexistent")
    assert result.returncode == 2
    assert "nonexistent" in result.stderr


# ---------------------------------------------------------------------------
# stats
# ---------------------------------------------------------------------------

def test_stats_exits_0_with_header(small_fastb):
    result = run("stats", small_fastb)
    assert result.returncode == 0
    assert "name\tlength\talpha\tenc\tgc_pct\tn_pct\tmasked_pct" in result.stdout
    assert "TOTAL" in result.stdout


# ---------------------------------------------------------------------------
# verify
# ---------------------------------------------------------------------------

def test_verify_good_file_exits_0(small_fastb):
    result = run("verify", small_fastb)
    assert result.returncode == 0
    assert "OK" in result.stdout


def test_verify_truncated_file_exits_3(small_fastb, tmp_path):
    bad = tmp_path / "truncated.fastb"
    with open(small_fastb, "rb") as f:
        bad.write_bytes(f.read(30))
    result = run("verify", str(bad))
    assert result.returncode == 3


# ---------------------------------------------------------------------------
# --version
# ---------------------------------------------------------------------------

def test_version():
    result = run("--version")
    assert "3.1.0" in result.stdout + result.stderr
