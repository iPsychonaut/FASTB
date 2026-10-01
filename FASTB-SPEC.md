# FASTB Specification, version 3.1

Status: draft. This file is the single authority. README and code follow it.
Date: 2026-10-01.

## 1. Purpose

FASTB stores DNA and RNA sequences at 2 bits per base. It replaces gzip on
FASTA intermediate files: smaller on disk, faster to read, no compression step.

The sequence core uses the same idea as the UCSC .2bit format: every base is
stored as 2 bits, and positions that hold anything else (N, ambiguity codes,
lowercase) are recorded as run lists next to the packed data. FASTB adds, on
top of .2bit:

| Addition | What it gives |
|---|---|
| Per-record CRC32 | Detects corrupted or truncated payloads |
| UTF-8 record headers of any length | No name length limit |
| RNA flag (`ALPHA=R`) | Round-trips U without a separate format |
| IUPAC exception list | Ambiguity codes and gaps without leaving 2-bit |
| Key=value header fields | Optional fields can be added later |
| Plain-text index footer | Jump to any record by name or number |

## 2. Design choices

Each choice: one line why, one line what it costs.

- 2-bit core for every base. Why: a 4-bit tier for the whole record doubles it whenever one base is not A, C, G, or T. Cost: every other symbol needs a side list.
- N and lowercase stored as run lists. Why: in assemblies these come in long runs, so a (start, length) pair costs a few bytes per run. Cost: a record with N or lowercase scattered one base at a time grows the header by about 6 bytes per run.
- IUPAC codes and gaps stored as an exception list. Why: polishers such as Pilon write a handful of ambiguity codes per contig; measured on a 4.66 Mb hybrid E. coli assembly, 6 such bases put 3 of 17 contigs in 4-bit and made the file 1,755,783 bytes against 1,167,801 with the list (pigz -6: 1,452,135). Cost: about 10 bytes of header per exception, so a record dense in ambiguity codes is still written at 4 bits (section 6.3).
- Text headers and index. Why: `grep '^>'`, `tail`, and `file` work without a parser. Cost: about 60 to 80 bytes per record more than a fixed binary header.
- Explicit `BYTES` count and a newline terminator. Why: a decoder never guesses where a record ends (the v2 CRC-probe heuristic could truncate records). Cost: none in practice.
- No file-level CRC by default (`FILE_CRC=00000000`). Why: streaming writers cannot know it up front. Cost: whole-file tampering is only caught per record.

## 3. File identification

- Extension: `.fastb`
- MIME type (proposed): `application/vnd.fastb`
- Magic: the first 12 bytes are exactly `##FASTB 3.1\n`. Readers of this
  version must also accept `##FASTB 3.0\n` (section 8).

## 4. File layout

```
##FASTB 3.1\n                 magic line
# ...\n                       zero or more file comment lines
<record 0>
<record 1>
...
<index entry per record>      NAME \t header_offset \t length \n
##INDEX\tENTRIES=<n>\tOFFSETS_START=<byte>\n
##END\tFILE_CRC=<hex8>\n
```

`header_offset` is the byte position of the `>` that starts that record's
header line. `length` is the record length in bases.

## 5. Record layout

```
# optional comment line(s), UTF-8, attached to the next record
>NAME\tLEN=<n>\tENC=<2|4>\tALPHA=<D|R>\tCRC=<hex8>\tBYTES=<n>[\tNRUNS=<runs>][\tIUPAC=<runs>][\tMASK=<runs>]\n
<BYTES bytes of payload>
\n
```

Fields, in this order:

| Field | Meaning |
|---|---|
| `NAME` | Record identifier. No tab, space, newline, or leading `>`. Other UTF-8 allowed. |
| `LEN` | Number of bases, decimal. |
| `ENC` | `2` or `4` bits per base. See section 6. |
| `ALPHA` | `D` for DNA (2-bit `11` decodes to T) or `R` for RNA (`11` decodes to U). `P` is reserved for protein; readers must reject it. |
| `CRC` | CRC32 (zlib polynomial) of the payload bytes, 8 lowercase hex digits. |
| `BYTES` | Exact payload byte count. Readers advance by this, never by scanning. |
| `NRUNS` | Optional. Positions that decode to N. Only written for `ENC=2`. |
| `IUPAC` | Optional. Positions that decode to an ambiguity code or gap. Only written for `ENC=2`. |
| `MASK` | Optional. Positions that decode as lowercase. |

`NRUNS` and `MASK` are lists of `start:length` in lowercase hex, 0-based,
comma separated: `NRUNS=10:a,26:8` means 10 N at position 16 and 8 N at
position 38.

`IUPAC` is a list of `start:length:symbol`, same hex, where `symbol` is one
uppercase letter from `W S M K R Y B D H V` or `-` for a gap:
`IUPAC=28:1:K,c8:2:R` means K at position 40 and RR at positions 200 and 201.

In all three lists runs are sorted by start, do not overlap, and adjacent
runs that could be one run (same symbol for `IUPAC`) are written as one.

Readers must ignore header keys they do not know. Writers must not reuse a
retired key name for a different meaning. A file may mix `ALPHA` values across
records.

## 6. Payload encoding

### 6.1 ENC=2

```
A = 00   C = 01   G = 10   T/U = 11
N and every IUPAC symbol = 00 (restored from NRUNS and IUPAC)
```

Four bases per byte, first base in the two most significant bits. The last
byte is zero padded in its low bits. `BYTES = ceil(LEN / 4)`.

Decoding order: unpack all bases, set every `NRUNS` range to N, set every
`IUPAC` range to its symbol, then lowercase every `MASK` range. So `n` and
`k` come back through `MASK`.

A and T, and C and G, are bitwise complements. Reverse complement is a byte
NOT plus a 2-bit reversal within each byte.

### 6.2 ENC=4

One-hot IUPAC nibble: bit 3 = A, bit 2 = C, bit 1 = G, bit 0 = T/U. A
degenerate code is the OR of the bases it stands for.

```
A=1000  C=0100  G=0010  T/U=0001
W=1001  S=0110  M=1100  K=0011  R=1010  Y=0101
B=0111  D=1011  H=1101  V=1110  N=1111  gap=0000
```

Two bases per byte, first in the high nibble. `BYTES = ceil(LEN / 2)`. Gap
decodes as `-`. `.` on input is written as gap in both encodings. `NRUNS` and
`IUPAC` are not written for ENC=4 records; readers apply them if present.

### 6.3 Choosing ENC

Writers must make the same choice so that every implementation produces the
same bytes. Build the `IUPAC` list text `T` for the record (empty if the
record holds only A, C, G, T/U, N). Then:

- If `T` is empty, write ENC=2.
- Else if `7 + len(T) <= ceil(LEN / 2) - ceil(LEN / 4)`, write ENC=2 with
  `IUPAC=T`. The 7 is the length of `\tIUPAC=`; the right side is the number
  of bytes 4-bit packing would add.
- Else write ENC=4.

Why: the rule picks whichever is smaller on disk. What it costs: a writer
must build the list before it knows the encoding. A writer may stop early
once the list is certain to exceed the bound.

## 7. Reader rules

1. Read the magic line. Accept `##FASTB 3.1` and `##FASTB 3.0`. Reject
   anything else, including higher version numbers.
2. Skip `#` lines before the first record. A `#` line between records is a
   comment for the next record.
3. On a `>` line: parse fields, read exactly `BYTES` bytes, check CRC, read
   one byte and require `\n`. Decode.
4. On a line that is neither `#` nor `>` nor blank: this is the first index
   entry. Stop iterating records.
5. `##INDEX` or `##END`: stop.

To jump to record k or to a named record: seek to end of file, read the last
4 KB, find the last `\n##INDEX\t`, parse `OFFSETS_START` and `ENTRIES`, seek
to `OFFSETS_START`, and read `ENTRIES` lines. Entry k is record k.

If the index is missing or damaged, records can still be read from the top.

## 8. Versioning

The version is the magic line. A change that would make an older reader
decode wrong bases requires a new version number, new fixtures, and a
changelog entry in one commit. An optional key that only adds information
does not.

| Version | What changed | Readable by a 3.1 reader |
|---|---|---|
| 3.0 | First v3 layout. Any IUPAC code other than N forced ENC=4. | yes |
| 3.1 | `IUPAC=` exception list; ENC chosen by the rule in 6.3. | yes |

A 3.0 reader refuses a 3.1 file at the magic line. That is intended: it
would otherwise ignore `IUPAC=` and print A where the ambiguity code belongs.

Versions 1 and 2 were exploratory formats that no released reader supports.

## 9. Conformance

`tests/fixtures/` holds FASTA inputs and their `.fastb` outputs. Every
implementation must produce those exact bytes from those inputs and decode
those files back to the same sequences. Byte equality is the test.
`tests/fixtures/v3.0/` holds 3.0 files that every reader must still decode.

## 10. Out of scope

Quality scores (FASTQ), per-position annotations, compression, encryption.
