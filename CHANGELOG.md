# Changelog

## [3.0.0] - UNRELEASED

Initial public release.

### Added
- FASTB v3 format: ASCII-header + bit-packed binary payload hybrid.
- O(1) random access by record name via footer index.
- CRC32 payload integrity per record.
- `fastb` CLI with subcommands: view, head, stats, extract, encode, decode, verify.
- Alphabet detection layer (`fastb/alphabet.py`) with explicit extension point for future amino-acid support.
- Reserved `ALPHA=AA` header value for future protein records.
- Specification document at [docs/FASTB_v3_spec.md](docs/FASTB_v3_spec.md).

### Notes
- Previous FASTB versions (v1, v2, v2.1) were exploratory and are not preserved. v3 is the first production release.
- `ALPHA=` header field replaces the previous `NUC=` field name, using full-word values (DNA, RNA) instead of single letters (D, R) for readability and forward compatibility with amino-acid records.
