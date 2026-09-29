<div align="center">
  <img src="resources/FASTB_banner.png" alt="FASTB Banner" width="1000">
</div>

# FASTB

A 2-bit binary format for DNA and RNA. Smaller than gzipped FASTA, faster to
read, and no compression step. Built to replace `pigz` on the FASTA
intermediate files that [EGAP](https://github.com/iPsychonaut/EGAP) writes.

Format definition: [FASTB-SPEC.md](FASTB-SPEC.md). Changes: [CHANGELOG.md](CHANGELOG.md).

## Install

```bash
pip install fastb
```

From a checkout, `pip install -e .[dev]` adds pytest. `pip install fastb[bench]`
adds psutil for the benchmark script. A Bioconda recipe is in `recipe/`.

Requires Python 3.9 or newer and numpy.

## Command line

| Command | What it does |
|---|---|
| `fastb encode in.fasta [-o out.fastb]` | FASTA to FASTB. Accepts `.fasta`, `.fa`, `.fna`, and any of those with `.gz`. |
| `fastb decode in.fastb [-o out.fasta]` | FASTB to FASTA, 80 columns (`-w 0` for one line per record). |
| `fastb cat in.fastb` | FASTA to stdout. For pipes and process substitution. |
| `fastb head in.fastb [-n 5]` | First N records as FASTA. |
| `fastb extract in.fastb NAME...` | Named records, read through the index without scanning the file. |
| `fastb index in.fastb` | Record number, name, byte offset, length. |
| `fastb info in.fastb` | Record count, bases, bytes, bits per base, 2-bit and 4-bit record counts. |
| `fastb stats in.fastb` | Per-record GC, N, and lowercase percentages. |
| `fastb verify in.fastb` | Check magic line, every CRC, and the index. |

Feeding a tool that reads plain FASTA:

```bash
mafft <(fastb cat genome.fastb) > aligned.fasta
```

Exit codes: 0 ok, 1 error, 2 some records missing, 3 bad file, 4 not nucleotide.

`fastb encode` rejects protein FASTA. Sequences with any of `F I L P Q E Z`
are rejected outright. Sequences of 20 or more bases with more than 5% IUPAC
degenerate codes (`W S M K R Y B D H V`) are also rejected; pass
`--force-nucleotide` if that is a real nucleotide record.

## Python

```python
import fastb

fastb.write_records("out.fastb", [("chr1", "ACGTNNNNacgt")])   # streams, never buffers the file
for name, seq in fastb.iter_records("out.fastb"):             # one record at a time
    ...
rec = fastb.read_record("out.fastb", "chr1")                   # by name, through the index
rec = fastb.read_record("out.fastb", 0)                        # by position
```

## How it stores sequence

Every base is 2 bits (A=00, C=01, G=10, T/U=11). Positions holding N, and
positions that were lowercase, are stored as (start, length) run lists in the
record header. This is the UCSC `.2bit` idea. A record only falls back to 4
bits per base when it holds an IUPAC code other than N, or a gap.

Each record carries its length, encoding, alphabet (DNA or RNA), CRC32, and
payload byte count in a plain-text header line, so `grep '^>'` works on a
`.fastb` file. A plain-text index at the end of the file gives direct access
to any record.

## Measurements

Synthetic 20 Mb single-record FASTA, 80 columns, 10% of bases in lowercase
runs and 1% in N runs. Windows 11, Python 3.14, numpy 2.5, one thread.
Command: `python bench/bench.py syn.fasta`. gzip is used because pigz is not
on this machine; both are single-threaded here.

| Tool | Bytes on disk | Encode s | Decode s | Peak RSS MB |
|---|---|---|---|---|
| gzip -6 | 6,460,810 | 2.97 | 0.11 | 9 |
| fastb | 5,022,303 | 0.46 | 0.28 | 138 |

Of the 0.28 s decode, 0.06 s is decoding and the rest is interpreter and
numpy start-up, so the gap closes on larger files. The Phase 4 gate (three
real EGAP intermediate files against `pigz -6 -p N` and `pigz -dc`) has not
been run yet; its table will replace this one.

## Not in scope

Quality scores (FASTQ), FAST5, per-position annotations, compression.

## License

MIT. See [LICENSE](LICENSE).
