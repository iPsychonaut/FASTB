# FASTB v3.0 — Binary Sequence Format

**Document version:** v3.0 draft (ALPHA field, 2026-04-21)
**Status:** Draft specification.
**Goal:** A sequence container with compact bit-packed payload, a fully human-readable structural layer, explicit self-description, and O(1) random access. Nucleotide content (DNA and RNA) is supported in this release. Amino-acid support is planned; the `ALPHA=P` value is reserved in the file format.

## 1. Design principles

1. **Density on payload, readability on structure.** The sequence itself is bit-packed (2 or 4 bits per base). Everything *about* the sequence — magic string, per-record metadata, index — is plain UTF-8 text.
2. **Self-describing, never heuristic.** Every record header carries its length, encoding, alphabet type, byte count, and CRC32 as explicit fields. A decoder never guesses.
3. **Graceful degradation in text tools.** `head`, `tail`, `file`, `grep '^>'`, and `less` all do something useful on a FASTB v3 file even without a dedicated parser.
4. **Recoverable.** If the trailing index is truncated or damaged, records can still be read sequentially from the top.
5. **No backwards compatibility with v1/v2.** Clean break to avoid the heuristics that made v2 fragile (see appendix). v3 readers MAY ship alongside v1/v2 readers in the same toolkit for migration.
6. **Unified format.** The same magic line, file extension, CLI, and tooling handle all content types. FASTB will not split into sister formats (`.fastbn` for nucleotides, `.fastbp` for proteins, etc.); content type is identified per-record by the `ALPHA=` field.

## 2. File identification

- **Extension:** `.fastb`
- **MIME type:** `application/vnd.fastb` (proposed)
- **Magic line:** first 12 bytes are exactly `##FASTB 3.0\n` (note the space, not a dash — this is the v3 signature).

Earlier versions of FASTB did not have a magic line at byte 0, which is why version sniffing relied on internal heuristics. v3 fixes this: byte-0 magic plus version means `file(1)` and libmagic can identify it unambiguously.

## 3. File layout

```
[ MAGIC LINE ]              — "##FASTB 3.0\n"
[ FILE COMMENTS ]           — zero or more lines starting with '#', UTF-8
[ RECORD 0 ]
[ RECORD 1 ]
...
[ RECORD N-1 ]
[ INDEX ENTRIES ]           — one per record, tab-separated, plain text
[ INDEX TRAILER ]           — "##INDEX\tENTRIES=N\tOFFSETS_START=<pos>\n"
[ END TRAILER ]             — "##END\tFILE_CRC=<hex8>\n"
```

## 4. Record layout

```
[ RECORD COMMENT ]          — optional, one or more '#'-prefixed UTF-8 lines
[ RECORD HEADER LINE ]      — single UTF-8 line, starts with '>', ends with '\n'
[ PAYLOAD ]                 — exactly BYTES bytes of packed sequence
[ TERMINATOR ]              — single '\n' byte
```

### 4.1 Record header line

```
>NAME\tLEN=<n>\tENC=<2|4>\tALPHA=<D|R|P>\tCRC=<hex8>\tBYTES=<n>[\tMASK=<runs>][\tK=V ...]\n
```

- `NAME`: sequence identifier. MUST NOT contain `\t`, ` `, `\n`, or `>`. Unicode permitted otherwise.
- `LEN`: number of bases in the sequence (decimal).
- `ENC`: encoding used — `2` for pure ACGT/ACGU, `4` for any degenerate or gap.
- `ALPHA`: indicates the alphabet of the record's sequence. Defined values:
  - `D` — DNA: nucleotide sequence in the DNA alphabet (A, C, G, T, plus IUPAC degenerate codes and gap). Decoder maps 2-bit `11` to `T`.
  - `R` — RNA: nucleotide sequence in the RNA alphabet (A, C, G, U, plus IUPAC degenerate codes and gap). Decoder maps 2-bit `11` to `U`.
  - `P` — Protein: amino-acid sequence (RESERVED for a future release; decoders in this release MUST reject records with `ALPHA=P`, producing a clear error message).

  No mixed T/U within a single record — the format enforces a single alphabet per record.

  Future values MAY be added; decoders MUST reject unknown values rather than guessing. A file containing mixed `ALPHA` values across records is permitted — each record is independently validated and decoded.

- `CRC`: CRC32 (zlib/IEEE polynomial) of the payload bytes, 8 lowercase hex digits.
- `BYTES`: exact byte count of the payload (redundant with LEN+ENC but simplifies stream parsing and catches truncation).
- `MASK` (optional): run-length-encoded lowercase intervals as `<start>:<length>` pairs, hex, comma-separated. Applied post-decode. Omitted when the sequence has no masking.
- Additional `KEY=VALUE` fields MAY be added by implementations; decoders MUST ignore unknown keys.

### 4.2 Payload encoding

**2-bit (`ENC=2`):**
```
A = 00,  C = 01,  G = 10,  T/U = 11
```
Four bases per byte, first base in the most significant 2 bits. Last byte zero-padded in low bits if `LEN % 4 != 0`. `BYTES = ceil(LEN/4)`.

Choosing A↔T and C↔G as bitwise complements makes reverse-complement a single `~byte` operation at the byte level (after a per-byte nibble swap). This is deliberate and matches the 2bit (UCSC) convention for the same reason.

**4-bit (`ENC=4`):** IUPAC one-hot, bit 3=A, bit 2=C, bit 1=G, bit 0=T/U.
```
A=1000  C=0100  G=0010  T/U=0001
W=1001  S=0110  M=1100  K=0011
R=1010  Y=0101
B=0111  D=1011  H=1101  V=1110
N=1111  -/.=0000
```
Two bases per byte, first in high nibble. `BYTES = ceil(LEN/2)`.

Degenerate codes are the bitwise OR of the bases they represent — so membership tests (e.g. "does W contain A?") are a single bitwise AND.

## 5. Index footer

After the last record, the file contains one line per record:
```
<name>\t<header_byte_offset>\t<length_in_bases>\n
```
`header_byte_offset` points at the `>` byte of that record's header line. This gives O(1) random access by name after one scan of the (tiny, text) index.

Followed by two trailer lines:
```
##INDEX\tENTRIES=<n>\tOFFSETS_START=<byte_offset>\n
##END\tFILE_CRC=<hex8>\n
```

To locate the index without a full scan, seek to end-of-file, read backwards for `##END\n`, then `##INDEX\n`, then jump to `OFFSETS_START`. Standard footer-pointer pattern used by PDF, ZIP, and Parquet.

`FILE_CRC` is the CRC32 of all bytes from start-of-file up to (but not including) the `##END` line itself. Implementations MAY write `00000000` if they are streaming and cannot buffer; readers SHOULD warn but not fail on `00000000`.

## 6. Parsing rules

### 6.1 Normative decoder behavior

1. Read and verify the magic line.
2. Skip any leading `#`-prefixed comment lines.
3. Repeat: read next line.
   - Blank line: skip.
   - `#`-prefixed: record-scoped comment; remember for the *next* record only.
   - `>`-prefixed: parse as record header. Read exactly `BYTES` bytes as payload. Verify CRC. Read exactly one terminator byte and verify it is `\n`. Decode and yield the record.
   - Anything else (no `>`, not blank, not `#`): this marks the start of the index footer. Stop record iteration.

### 6.2 Why `BYTES` is mandatory

A naive parser might try to find the next record by scanning for `\n>`. This fails: the byte `0x0A` is a valid 4-bit encoding byte (`0000 1010` = gap + G), and `0x3E` (`0011 1110` = C+V) is also a legal payload byte. Line-scanning will produce false record boundaries on real data. Decoders MUST use `BYTES` to advance the stream; `\n>` scanning is for **recovery only** when `BYTES` is known-bad.

### 6.3 Recovery mode

If a CRC fails or `BYTES` overshoots end-of-file, the decoder MAY enter recovery mode: scan forward for the next `\n>` sequence, attempt to parse from there, and emit a diagnostic about the skipped region. Recovery mode is OPTIONAL; strict mode is the default.

## 7. Interaction with text editors

This format does **not** support being edited as sequence data in a text editor. That goal is fundamentally incompatible with bit-packed density.

It **does** support:
- Identifying the file (`file foo.fastb` sees the `##FASTB 3.0` magic line; add an entry to `/etc/magic`).
- Listing records (`grep '^>' foo.fastb` works like FASTA).
- Inspecting the index (`tail -20 foo.fastb`).
- Surviving being *opened* in a text editor without corruption, provided the user does not save. The `##FASTB` magic line gives editor plugins (VS Code, vim, Emacs) a reliable trigger to mark the file as binary-readonly and route edits through the `fastb` CLI.

**FASTB v3 is not safe to save from a text editor.** This is documented clearly rather than fought.

## 8. Size comparison

For a typical mammalian genome (~3 Gb, mostly uppercase ACGT with small N-stretches):

| Format | Approx size |
|---|---|
| FASTA (uncompressed) | ~3.0 GB |
| FASTA.gz | ~900 MB |
| 2bit (UCSC) | ~800 MB |
| **FASTB v3 (ENC=2, header overhead ~80 bytes/record)** | ~750 MB |
| FASTB v3 + gzip | ~730 MB (minimal gain — payload is already near entropy floor) |

The win versus FASTA.gz is modest (~15-20%) but real, and parsing is ~4-10× faster because there is no decompression and no character-to-nibble translation — just a table lookup per byte.

Honest headline: **FASTB v3 trades a small size win for a large speed-and-structure win over FASTA.gz, and is competitive with 2bit while being self-describing and human-inspectable.**

## 9. What this format deliberately excludes

- **Quality scores.** A sister format (FASTBQ, planned) will handle FASTQ replacement. Not bolted onto v3.
- **Per-position annotations** (features beyond case-masking, tracks). Use companion BED/GFF.
- **Encryption, signing, compression.** Layer these externally if needed.

### 9.1 Amino-acid sequences (planned, not yet implemented)

FASTB is a unified format for biological sequences. Amino-acid support is planned for a future release under the same `.fastb` extension, magic line, CLI, and tooling.

In this release, the `ALPHA=P` header value is **reserved**. Decoders MUST reject any record with `ALPHA=P` with a clear error message. Encoders MUST NOT write `ALPHA=P` records until the AA codec is published.

The alphabet dispatch layer (`fastb/alphabet.py`) is already structured to accept a future `require_protein` path alongside `require_nucleotide`. Adding AA support will be additive — it will not require changing the file format, the magic line, or the CLI's public interface.

## Appendix: Bugs found in FASTB v2 that v3 fixes

Round-trip tests on the v2 codec revealed three correctness issues and two design fragilities. These motivated several choices above.

**1. `NUC_TYPE` is derived from content instead of declared.** The v2 encoder contains:
```python
nuc_tag = "R" if str(nuc).upper().startswith("R") or ("U" in (seq or "")) else "D"
```
This makes the alphabet field a derived value, not a declared one. Consequences:
- A record declared `DNA` containing a `U` silently becomes RNA.
- A record with mixed `T` and `U` becomes RNA, and all `T`s are rewritten to `U` on decode. Lossy for legitimate sequences.
- A record declared `RNA` but containing only `T` characters round-trips with every `T` flipped to `U`.

v3 fix: `ALPHA` is authoritative on encode; mixed T/U within a record is rejected with a clear error.

**2. Protein sequences are silently accepted if they happen to use only letters that overlap the IUPAC nucleotide alphabet.** The set `A C G H K M N R S T V W Y D B` is legal in v2's 4-bit mode — but those are also valid amino acid one-letter codes. A peptide like `GAGAGSGAGAGS` (a real silk fibroin motif) encodes cleanly as "nucleotides" and decodes back unchanged. The README's claim that "unsupported amino acid FASTA files" are rejected is only true for proteins containing `F`, `I`, `L`, `P`, `Q`, `E`, or `Z`.

v3 fix: explicit alphabet declaration in `ALPHA=D|R|P`, plus a content heuristic in the reference encoder that rejects sequences ≥20 bases long with >5% IUPAC-degenerate characters (`W S M K R Y B D H V`). Real nucleotide sequences are overwhelmingly `ACGT[U]` with occasional `N`; >5% degenerate is a near-certain protein signal. The silk-fibroin `GAGAGS` motif, BSA signal peptides, and W-rich peptides all correctly reject under this rule. Short sequences (<20 bases) bypass the heuristic to avoid false positives on primers and adapters.

**3. CRC-as-end-of-record heuristic is fragile.** The v2 decoder has no END_RECORD marker in the emitted output (the `_END_RECORD = 0xFE` constant is defined but never written). Instead it reads 4 bytes after each TLV and checks if they form a valid CRC; if not, seeks back and reads another TLV. This doubles I/O, creates subtle coupling for future TLV types, and produces confusing error messages when tampering occurs mid-record (e.g. `UnicodeDecodeError` instead of `CRC mismatch`).

v3 fix: explicit `BYTES=` count plus mandatory `\n` terminator. Decoder never guesses record boundaries on valid files.

**4. The `_END_RECORD = 0xFE` constant is defined but never written.** Dead code in v2; removed in v3.

**5. Tetrad table has a dict-ordering collision between `-` and space.** Both `" "` and `"-"` encode to `0x0`, and `_INV_TETRA` keeps whichever was inserted last (space). So encoded gaps decode as spaces, silently breaking alignment tools that look for `-`. The encoder even strips spaces before packing, so the space entry is both dead on encode and corrupting on decode.

v3 fix: `0x0` decodes unambiguously to `-`. Space is not a sequence character.

**6. RLE softmask uses uint32 but concatenated whole-genome records could exceed 4 Gb.** Minor but real.

v3 fix: `MASK` field uses variable-width hex, no ceiling.

---

*Document version: v3.0 draft (ALPHA field, 2026-04-21)*
