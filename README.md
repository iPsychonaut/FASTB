<div align="center">
  <img src="resources/FASTB_banner.png" alt="FASTB Banner" width="1000">
</div>

# FASTB

Binary sequence format. Compact, self-describing, and fast. Today for DNA and RNA; amino-acid support is planned under the same format.

---

## What it is

FASTB v3 is a sequence container with bit-packed payloads and a human-readable structural layer. Each record has an ASCII header (readable with `grep`, `tail`, or any text editor) and a compact binary payload. A footer index provides O(1) random access by record name.

---

## Why not FASTA

| | FASTA.gz | FASTB v3 |
|---|---|---|
| Size (typical mammalian genome) | ~900 MB | ~750 MB |
| Parsing | Decompress + ASCII→code | Table lookup per byte |
| Speed | Baseline | ~4–10× faster |
| Random access by name | Requires index tool | Built-in footer index |
| Self-describing | No | Yes (ASCII header per record) |

The size win is real but modest. The speed and structure win is substantial.

---

## Quick start

```bash
pip install -e .

fastb encode genome.fasta          # → genome.fastb
fastb head genome.fastb            # first 5 records as FASTA
fastb extract genome.fastb chr1    # random-access by name
fastb stats genome.fastb           # per-record statistics
fastb verify genome.fastb          # integrity check
fastb decode genome.fastb          # back to FASTA
```

---

## What's in the file

```
##FASTB 3.0                              ← magic line, identifies format
# optional comments
>chr1    LEN=248956422  ENC=2  ALPHA=D  CRC=4a3c1b2d  BYTES=62239106
<62239106 payload bytes>
>chr2    ...
...
chr1    85    248956422                  ← footer index (plain text)
chr2    ...
##INDEX    ENTRIES=24  OFFSETS_START=...
##END      FILE_CRC=...
```

Full specification: [docs/FASTB_v3_spec.md](docs/FASTB_v3_spec.md)

---

## CLI subcommands

| Subcommand | Description |
|---|---|
| `encode <file.fasta>` | Convert FASTA → FASTB. Detects and rejects protein sequences. |
| `decode <file.fastb>` | Convert FASTB → FASTA. |
| `view <file.fastb>` | Stream records as FASTA. Supports `--names` filter and index-based access. |
| `head <file.fastb>` | Print the first N records (default 5). |
| `extract <file.fastb> <name>...` | Random-access extraction by record name. |
| `stats <file.fastb>` | Per-record statistics: length, alpha, GC%, N%, masked%. |
| `verify <file.fastb>` | Check magic, CRCs, and index consistency. |

Run `fastb <subcommand> --help` for options.

---

## Supported content

- **DNA**: A, C, G, T plus IUPAC degenerate codes (W S M K R Y B D H V N) and gap (`-`).
- **RNA**: A, C, G, U plus IUPAC degenerate codes and gap.
- **Amino acids** (`ALPHA=P`): detected and rejected with a clear error (exit code 4). The file format reserves `ALPHA=P` for protein records; protein codec support is planned for a future release under the same `.fastb` extension, CLI, and tooling. FASTB is a unified format — there will be no sister protein format.

---

## Alphabet detection

`fastb encode` classifies each record's alphabet before writing:

- **Hard reject**: sequences containing `F I L P Q E Z` are amino acids (these letters are not in the IUPAC nucleotide alphabet). Exit code 4.
- **Heuristic reject**: sequences ≥20 bases long with >5% IUPAC-degenerate codes (`W S M K R Y B D H V`, excluding N) are almost certainly proteins (e.g. silk-fibroin motifs). Exit code 4. Bypassable with `--force-nucleotide`.
- `--force-nucleotide` skips the heuristic but never bypasses the hard-letter check.

---

## Roadmap

- Amino-acid sequence support (ALPHA=P codec, same format and CLI)
- FASTQ equivalent with quality scores
- Nanopore FAST5 integration

---

## License

MIT — see [LICENSE](LICENSE).
