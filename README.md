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

Requires Python 3.8 or newer and numpy.

## Command line

| Command | What it does |
|---|---|
| `fastb encode in.fasta [-o out.fastb] [--verify]` | FASTA to FASTB. Accepts `.fasta`, `.fa`, `.fna`, and any of those with `.gz`. `--verify` re-reads the output and compares every header line and sequence with the input; on a difference it removes the output and exits 3. |
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

Every base is 2 bits (A=00, C=01, G=10, T/U=11). Positions holding N,
positions holding an ambiguity code (K, R, Y and the other IUPAC letters) or
a gap, and positions that were lowercase are stored as short lists in the
record header. This is the UCSC `.2bit` idea, extended to ambiguity codes. A
record is written at 4 bits per base only when its ambiguity list would take
more room than 4-bit packing, which means roughly more than 1 base in 40.

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

## Inside EGAP

Three E. coli samples (ONT + Illumina hybrid, Illumina with a reference,
PacBio) run end to end through [EGAP](https://github.com/iPsychonaut/EGAP)
with pigz and with FASTB, twice: once in WSL and once on a local Kubernetes
cluster.

### Cluster run, format 3.1

k3d cluster on the machine above, one pod per sample at 7 CPUs and 24Gi, the
pigz and FASTB variants side by side. Images `egap:old` and `egap:fastb`
(EGAP at 543be05, FASTB 3.1.0); the FASTB Job sets
`EGAP_INTERMEDIATE_FORMAT=fastb`. Each pod undoes any compression in the
sample folder, then times `final_compress.py` alone and prints one
`COMPRESS_SUMMARY` line (EGAP `tests/ab/ab-test.yaml`). File sizes from
`find -printf '%s'`. Both Jobs completed 3 of 3 samples (3 h 37 min pigz,
3 h 44 min FASTB) with every assembler running, Flye included.

| Sample | Measure | pigz | FASTB 3.1 |
|---|---|---|---|
| Illumina | Final assembly | 1,443,530 B | 1,160,593 B |
| | All compressed FASTA | 63,193,972 B | 58,200,090 B |
| | `final_compress` time | 12.1 s | 29.0 s |
| Hybrid | Final assembly | 1,451,758 B | 1,168,250 B |
| | All compressed FASTA | 44,132,949 B | 40,264,419 B |
| | `final_compress` time | 14.3 s | 20.4 s |
| PacBio | Final assembly | 1,695,559 B | 1,364,249 B |
| | All compressed FASTA | 49,647,840 B | 44,498,569 B |
| | `final_compress` time | 6.2 s | 15.7 s |
| All three | Compressed FASTA | 156,974,761 B | 142,963,078 B |
| | `final_compress` time | 32.6 s | 65.2 s |

FASTB wrote 57 files, all at 2 bits per base; 8 carry an IUPAC list. The
assemblies have the same size, contig count, and N50 in both variants
(Illumina 4.64 Mb in 1 contig; hybrid 4.66 Mb, 17 contigs, N50 465 kb;
PacBio 5.46 Mb, 3 contigs, N50 2.94 Mb). Sample folders after compression
differ by up to 17 MB between variants, all of it in MaSuRCA's working
files, which vary run to run.

### WSL run, format 3.0 then 3.1

Same machine, 16 threads, 48 GB, `--intermediate_format pigz` then `fastb`,
each arm with its own fresh input and output folders, disk use sampled every
30 s with `du -sb`. The FASTB arm ran format 3.0; the 3.1 column re-encodes
that arm's 54 `.fastb` files afterwards (`fastb cat -w 0 x.fastb | fastb
encode`), with every round trip checked by `cmp`.

| | pigz | fastb 3.0 | fastb 3.1 |
|---|---|---|---|
| Wall time, 3 samples | 6,852 s | 6,843 s | not rerun |
| Peak disk | 23.523 GB | 23.517 GB | not rerun |
| Compressed FASTA, all files | 138,976,865 B | 129,496,900 B | 125,362,296 B |
| The 54 files FASTB encodes | 79,692,469 B | 70,243,637 B | 66,109,033 B |
| Hybrid final assembly (4.66 Mb, 6 IUPAC bases) | 1,452,135 B | 1,755,783 B | 1,167,801 B |

### What the numbers say

- Each final assembly is about 20% smaller as FASTB 3.1 than as pigz -6
  (19.5 to 19.6% in all three cluster samples).
- Across all of a run's FASTA files the saving is about 9%, because about
  half of them stay on pigz (see below). On a 5 Mb genome that is 14 MB out
  of several GB: reads, VCF, and GFF files are most of what EGAP writes.
- FASTB's compression step takes about twice as long as pigz here: 65 s
  against 33 s over three samples. Each file is decoded and compared with
  the original before the original is deleted, and Python starts twice per
  file. Total run time did not change measurably (6,843 s against 6,852 s
  in WSL; 3 h 44 min against 3 h 37 min in the cluster, where pods share
  the host).
- On the files it encodes, FASTB 3.1 is 17% smaller than pigz -6.
- Format 3.0 was larger than pigz on every Pilon-polished assembly: 6
  ambiguity codes put 3 of 17 contigs in 4-bit. That is why 3.1 has the
  IUPAC list.
- 46 of 100 FASTA files stayed on pigz, all for sound reasons: 36 protein
  files from Compleasm, and 10 that were empty or whose headers carry
  descriptions that FASTB would drop (Racon's `LN:i:` tags, for example).
- The final assemblies were not byte-identical between runs, with or without
  FASTB. The hybrid assembly differed between two pigz runs, and the one
  PacBio difference (4 bp in one contig) was already present in Flye's raw
  output, before FASTB is first used.

## Not in scope

Quality scores (FASTQ), FAST5, per-position annotations, compression.

## License

MIT. See [LICENSE](LICENSE).
