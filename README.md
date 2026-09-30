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

Synthetic single-record FASTA, 80 columns, 10% of bases in lowercase runs
and 1% in N runs. WSL2 Ubuntu 24.04 on a Ryzen 7 5800H (16 threads), pigz
2.8, Python 3.12, numpy 2.5.3, files on ext4, best of 3. Wall time and peak
RSS from `/usr/bin/time -f "%e %M"`. `fastb cat` output was byte-identical
to `pigz -dc` for every file. Script: `bench/sweep.sh`.

| Input | Tool | Bytes on disk | Encode s | Decode s | Peak RSS MB |
|---|---|---|---|---|---|
| 5 MB | pigz -6 -p 16 | 1,613,641 | 0.06 | 0.02 | 10 |
| 5 MB | fastb | 1,255,465 | 0.20 | 0.18 | 89 |
| 50 MB | pigz -6 -p 16 | 16,141,509 | 0.51 | 0.17 | 10 |
| 50 MB | fastb | 12,558,640 | 0.71 | 0.36 | 274 |
| 200 MB | pigz -6 -p 1 | 64,554,815 | 16.97 | 0.77 | 3 |
| 200 MB | pigz -6 -p 16 | 64,554,815 | 1.85 | 0.63 | 10 |
| 200 MB | fastb | 50,239,522 | 2.14 | 0.90 | 462 |
| 200 MB | fastb -p 4 | 50,239,522 | 2.14 | 0.75 | 121 |

What the numbers say:

- Bytes: fastb is 22% smaller than pigz -6 at every size.
- `pigz -dc` does not get faster with more threads (0.77 s at 1, 0.63 s at
  16). Inflate is serial. So the decode target is fixed at about 300 MB/s
  of FASTA out.
- fastb decode is about 0.2 s behind pigz at every size. That 0.2 s is
  Python and numpy start-up, not decoding: the 2-bit unpack of 200 MB takes
  0.11 s, wrapping 0.13 s, and writing the file 0.3 s (a plain `cat` of the
  same file takes 0.32 s on this machine). On files under about 10 MB the
  start-up cost dominates and pigz wins outright.
- `fastb -p N` splits records into 16 Mb chunks across N processes. It helps
  on large files (0.90 to 0.75 s) and hurts on small ones (process start).
- Encode: fastb single-threaded beats pigz up to 4 threads and loses to
  pigz -p 16 by about 15%.

The Phase 4 gate (three real EGAP intermediate files) has not been run yet.

## Not in scope

Quality scores (FASTQ), FAST5, per-position annotations, compression.

## License

MIT. See [LICENSE](LICENSE).
