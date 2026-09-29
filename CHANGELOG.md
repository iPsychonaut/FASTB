# Changelog

## [3.0.0] - UNRELEASED

First release. Format defined in [FASTB-SPEC.md](FASTB-SPEC.md).

### Format
- Every base is stored at 2 bits. N and lowercase positions are run lists in
  the header (`NRUNS=`, `MASK=`). 4-bit storage is used only for records with
  IUPAC codes other than N, or gaps.
- Records end with an explicit `BYTES=` count and a newline. No guessing.
- Plain-text index footer. Jump to record k or to a record by name.
- Readers ignore header keys they do not know.

### Code
- Encoder and decoder run on numpy lookup tables.
- `fastb` CLI: view, head, stats, extract, encode, decode, verify.
- Alphabet checks reject protein FASTA with exit code 4.
- Byte-exact fixtures in `tests/fixtures/`.

### Removed
- The v1 reader. It never returned data.
- The v2 TLV format, its CRC-probe end-of-record heuristic, and the 3-bit
  tier. No released reader supports v1 or v2 files.
- The HTML editor and the PyQt5 compare tool.
