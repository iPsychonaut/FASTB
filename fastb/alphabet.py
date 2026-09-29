"""Alphabet detection and validation for FASTB.

FASTB v3 currently encodes nucleotide content only. This module is
the classification and validation layer that will eventually dispatch
between nucleotide and amino-acid codecs. For now, `require_nucleotide`
is the only callable path to encoding; amino-acid sequences raise
ProteinDetectedError with a clear message.

Planned (not implemented in this release):
  - require_protein(): validation counterpart for AA sequences.
  - classify_and_route(): dispatcher that selects codec+validator
    based on detect_alphabet() output.

Design intent: adding AA support is additive. It does NOT require
changing detect_alphabet's signature, require_nucleotide's behavior,
or the public API of fastb/core.py. The ALPHA= header field was
designed for this extension.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Tunable thresholds
# ---------------------------------------------------------------------------

PROTEIN_HEURISTIC_MIN_LENGTH = 20
PROTEIN_HEURISTIC_DEGENERATE_FRACTION = 0.05

# Letters that are in the amino-acid alphabet but NOT in the IUPAC
# nucleotide alphabet. Their presence is definitive.
_PROTEIN_HARD_LETTERS = frozenset("FILPQEZ")

# IUPAC degenerate codes that overlap with amino-acid one-letter codes.
# N is intentionally excluded: real genomic data is routinely 10%+ N.
_DEGENERATE_LETTERS = frozenset("WSMKRYBDHV")

_LEGAL_DNA = frozenset("ACGTWSMKRYBDHVN-.")
_LEGAL_RNA = frozenset("ACGUWSMKRYBDHVN-.")


# ---------------------------------------------------------------------------
# Exception hierarchy
# ---------------------------------------------------------------------------

class AlphabetError(ValueError):
    """Base class for alphabet-related errors. Exit code 4 from the CLI."""


class ProteinDetectedError(AlphabetError):
    """Sequence appears to be amino acids.

    Distinct from other AlphabetErrors so that scripts can route proteins
    to a future AA pipeline without parsing error messages.
    """


# ---------------------------------------------------------------------------
# Classifier
# ---------------------------------------------------------------------------

def detect_alphabet(seq_upper: str) -> str:
    """Classify a sequence as 'dna', 'rna', 'protein', or 'ambiguous'.

    Pure classifier. Does NOT raise. Callers decide what to do with
    the result.

    Returned strings are lowercase category labels used internally.
    They are NOT the same as the file-format ALPHA= values (which are
    single letters: D, R, P). This case split is intentional. It keeps
    "classifier decided" distinct from "header field value".

    Rules, in priority order:

      1. Contains any of F I L P Q E Z        -> 'protein'
         These letters are in the amino-acid alphabet but NOT in the
         IUPAC nucleotide alphabet.

      2. Contains both T and U                -> 'ambiguous'
         Nucleotide data uses T (DNA) or U (RNA) exclusively.

      3. len(seq) >= PROTEIN_HEURISTIC_MIN_LENGTH
           AND degenerate_fraction > PROTEIN_HEURISTIC_DEGENERATE_FRACTION
                                              -> 'protein'
         Catches IUPAC-overlap proteins (e.g. silk fibroin GAGAGS motif).
         Degenerate letters counted: W S M K R Y B D H V. N is excluded.

      4. U present (rule 2 didn't fire)       -> 'rna'
      5. T present                            -> 'dna'
      6. Otherwise                            -> 'ambiguous'
    """
    # Rule 1: definitive protein letters
    for ch in seq_upper:
        if ch in _PROTEIN_HARD_LETTERS:
            return "protein"

    has_t = "T" in seq_upper
    has_u = "U" in seq_upper

    # Rule 2: mixed T and U
    if has_t and has_u:
        return "ambiguous"

    # Rule 3: degenerate-fraction heuristic
    if len(seq_upper) >= PROTEIN_HEURISTIC_MIN_LENGTH:
        deg_count = sum(1 for ch in seq_upper if ch in _DEGENERATE_LETTERS)
        if deg_count / len(seq_upper) > PROTEIN_HEURISTIC_DEGENERATE_FRACTION:
            return "protein"

    # Rules 4-6
    if has_u:
        return "rna"
    if has_t:
        return "dna"
    return "ambiguous"


# ---------------------------------------------------------------------------
# Validation chokepoint
# ---------------------------------------------------------------------------

def _letter_counts(seq: str) -> "np.ndarray":
    """Count of each byte value, case-folded to uppercase. Shape (256,)."""
    raw = np.frombuffer(seq.encode("ascii", errors="replace"), dtype=np.uint8)
    # bincount upcasts its input to int64; chunking keeps that temp at 8 MB.
    counts = np.zeros(256, dtype=np.int64)
    for i in range(0, len(raw), 1 << 20):
        counts += np.bincount(raw[i:i + (1 << 20)], minlength=256)
    counts[ord("A"):ord("Z") + 1] += counts[ord("a"):ord("z") + 1]
    counts[ord("a"):ord("z") + 1] = 0
    return counts


def require_nucleotide(
    seq: str,
    declared_alpha: str,
    name: str,
    force_nucleotide: bool = False,
) -> None:
    """Validate that a sequence is encodable as nucleotide content.

    Raises ProteinDetectedError if the sequence classifies as protein.
    Raises AlphabetError for other mismatches: mixed T/U, declared
    ALPHA=DNA but sequence contains U, declared ALPHA=RNA but sequence
    contains T, letters outside the IUPAC nucleotide alphabet.

    This is the ONLY place in the core encoder that gates what gets
    written to a v3 file as nucleotide content. When amino-acid support
    lands, a sibling `require_protein` function will be added.

    Args:
        seq: sequence string, any case.
        declared_alpha: "D" or "R" from the Record.
        name: record name, used in error messages.
        force_nucleotide: if True, skip the degenerate-fraction heuristic
            (rule 3 in detect_alphabet). Does NOT bypass the hard-letter
            check (rule 1). F/I/L/P/Q/E/Z are never nucleotides.
    """
    counts = _letter_counts(seq)

    def present(letters) -> str:
        return "".join(ch for ch in sorted(letters) if counts[ord(ch)])

    # Hard protein letters: never nucleotides, no bypass possible
    hard = present(_PROTEIN_HARD_LETTERS)
    if hard:
        raise ProteinDetectedError(
            f"Record {name!r}: sequence appears to be amino acids "
            f"(contains {hard[0]!r} which is not in the IUPAC nucleotide alphabet). "
            f"FASTB v3 encodes nucleotides only; amino-acid support is planned "
            f"for a future release."
        )

    has_t = bool(counts[ord("T")])
    has_u = bool(counts[ord("U")])
    if has_t and has_u:
        raise AlphabetError(
            f"Record {name!r}: contains both T and U. FASTB v3 requires a single "
            f"nucleotide alphabet per record (ALPHA=D uses T, ALPHA=R uses U)."
        )
    if declared_alpha == "D" and has_u:
        raise AlphabetError(
            f"Record {name!r}: declared ALPHA=D but sequence contains U. "
            f"Declare ALPHA=R or convert U to T before encoding."
        )
    if declared_alpha == "R" and has_t:
        raise AlphabetError(
            f"Record {name!r}: declared ALPHA=R but sequence contains T. "
            f"Declare ALPHA=D or convert T to U before encoding."
        )

    # Illegal characters
    legal = _LEGAL_DNA if declared_alpha == "D" else _LEGAL_RNA
    seen = np.flatnonzero(counts)
    illegal = [chr(b) for b in seen if chr(b) not in legal]
    if illegal:
        raise AlphabetError(
            f"Record {name!r}: illegal symbol {illegal[0]!r} for ALPHA={declared_alpha}. "
            f"Legal characters: {sorted(legal)}"
        )

    # Degenerate-fraction heuristic (bypassable with --force-nucleotide)
    if not force_nucleotide and len(seq) >= PROTEIN_HEURISTIC_MIN_LENGTH:
        deg_count = int(sum(counts[ord(ch)] for ch in _DEGENERATE_LETTERS))
        fraction = deg_count / len(seq)
        if fraction > PROTEIN_HEURISTIC_DEGENERATE_FRACTION:
            raise ProteinDetectedError(
                f"Record {name!r}: sequence appears to be amino acids "
                f"({fraction:.1%} of residues are IUPAC-degenerate codes, above the "
                f"{PROTEIN_HEURISTIC_DEGENERATE_FRACTION:.1%} heuristic threshold). "
                f"FASTB v3 encodes nucleotides only; amino-acid support is planned "
                f"for a future release. If this is genuinely a heavily-degenerate "
                f"nucleotide sequence, re-run with --force-nucleotide."
            )
