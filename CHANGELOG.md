# Changelog

## [3.1.0] - 2026-10-07

First release. Writes file format 3.1, reads 3.1 and 3.0. Format defined in
[FASTB-SPEC.md](FASTB-SPEC.md).

### Format
- Every base is stored at 2 bits. N and lowercase positions are run lists in
  the header (`NRUNS=`, `MASK=`).
- 3.1: IUPAC ambiguity codes and gaps are an exception list (`IUPAC=`) over
  the 2-bit core. 4-bit storage is used only when that list would be larger
  than the bytes 4-bit packing adds. In 3.0 a single ambiguity code put the
  whole record in 4-bit; on a Pilon-polished 4.66 Mb assembly that was
  1,755,783 bytes against 1,167,801 now.
- The magic line is `##FASTB 3.1`. A 3.0 reader refuses a 3.1 file instead of
  misreading it. 3.0 files stay readable; fixtures for both are in
  `tests/fixtures/`.
- Records end with an explicit `BYTES=` count and a newline. No guessing.
- Plain-text index footer. Jump to record k or to a record by name.
- Readers ignore header keys they do not know.

### Code
- Encoder and decoder run on numpy lookup tables in 16 Mb chunks; memory
  is bounded per chunk, not per record.
- `fastb cat -p N` and `fastb decode -p N` decode chunks in N processes.
- `fastb` CLI: view, head, stats, extract, encode, decode, verify.
- `fastb encode --verify` re-reads the written file and compares every
  header line and sequence with the input; on a difference it removes the
  output and exits 3. A dropped header description counts as a difference.
- `fastb encode` takes several files in one process, each encoded and
  verified on its own; `--append` names outputs `<input>.fastb`; the exit
  code is the highest per-file code.
- Alphabet checks reject protein FASTA with exit code 4.
- Byte-exact fixtures in `tests/fixtures/`.

### Removed
- The v1 reader. It never returned data.
- The v2 TLV format, its CRC-probe end-of-record heuristic, and the 3-bit
  tier. No released reader supports v1 or v2 files.
- The HTML editor and the PyQt5 compare tool.
