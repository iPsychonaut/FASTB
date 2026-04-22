"""Tests for fastb.alphabet: classifier and validation chokepoint."""

import pytest

from fastb.alphabet import (
    PROTEIN_HEURISTIC_DEGENERATE_FRACTION,
    PROTEIN_HEURISTIC_MIN_LENGTH,
    AlphabetError,
    ProteinDetectedError,
    detect_alphabet,
    require_nucleotide,
)


# ---------------------------------------------------------------------------
# detect_alphabet classification table
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seq,expected", [
    ("ACGTACGTACGTACGTACGT",         "dna"),       # 20 pure ACGT
    ("ACGUACGUACGUACGUACGU",         "rna"),       # 20 pure ACGU
    ("GAGAGSGAGAGSGAGAGSGAGAGS",     "protein"),   # heuristic: 4/24 S > 5%
    ("MKWVTFISLLLLFSSAYS",           "protein"),   # hard letters F/I/L
    ("GAGAGS",                        "ambiguous"), # len=6 < 20 AND no T/U → rule 6
    # 'ACGT' has T → rule 5 fires ('dna'); rules 4-6 have no length floor.
    ("ACGT",                          "dna"),
    ("ACGTACGTACGTACGTACGTN",        "dna"),       # N doesn't count; rule 5 fires
    ("ACGTUACGT",                    "ambiguous"), # mixed T/U → rule 2
])
def test_detect_alphabet(seq, expected):
    assert detect_alphabet(seq) == expected


# ---------------------------------------------------------------------------
# require_nucleotide — happy path
# ---------------------------------------------------------------------------

def test_require_nucleotide_valid_dna():
    require_nucleotide("ACGTACGT", "D", "rec")  # should not raise


def test_require_nucleotide_valid_rna():
    require_nucleotide("ACGUACGU", "R", "rec")


def test_require_nucleotide_with_degenerates():
    require_nucleotide("ACGTWSMKRYBDHVN", "D", "rec")


# ---------------------------------------------------------------------------
# require_nucleotide — protein detection
# ---------------------------------------------------------------------------

def test_require_nucleotide_rejects_hard_letter():
    with pytest.raises(ProteinDetectedError, match="amino acid"):
        require_nucleotide("FILPQEZACGT", "D", "test_rec")


def test_require_nucleotide_rejects_heuristic_protein():
    # GAGAGSGAGAGS... repeated to be >= PROTEIN_HEURISTIC_MIN_LENGTH
    seq = "GAGAGS" * 5  # len=30; S fraction = 5/30 ≈ 16.7% > 5%
    with pytest.raises(ProteinDetectedError, match="heuristic"):
        require_nucleotide(seq, "D", "silk")


def test_force_nucleotide_bypasses_heuristic():
    """--force-nucleotide must skip the heuristic but NOT the hard-letter check."""
    seq = "GAGAGS" * 5
    require_nucleotide(seq, "D", "silk", force_nucleotide=True)  # should not raise

    # Hard letters still rejected even with force_nucleotide
    with pytest.raises(ProteinDetectedError):
        require_nucleotide("FILPQEZ" + seq, "D", "silk", force_nucleotide=True)


# ---------------------------------------------------------------------------
# require_nucleotide — alphabet mismatch
# ---------------------------------------------------------------------------

def test_require_nucleotide_rejects_mixed_tu():
    with pytest.raises(AlphabetError, match="both T and U"):
        require_nucleotide("ACGTUACGT", "D", "rec")


def test_require_nucleotide_rejects_dna_with_u():
    with pytest.raises(AlphabetError, match="ALPHA=D"):
        require_nucleotide("ACGUACGU", "D", "rec")


def test_require_nucleotide_rejects_rna_with_t():
    with pytest.raises(AlphabetError, match="ALPHA=R"):
        require_nucleotide("ACGTACGT", "R", "rec")


# ---------------------------------------------------------------------------
# Constants are importable at module scope
# ---------------------------------------------------------------------------

def test_constants_accessible():
    assert PROTEIN_HEURISTIC_MIN_LENGTH == 20
    assert PROTEIN_HEURISTIC_DEGENERATE_FRACTION == 0.05
