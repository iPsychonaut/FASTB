"""FASTB command-line interface.

Exit codes: 0 success, 1 generic error, 2 partial success (some records
missing), 3 file format error, 4 alphabet error.
"""

from __future__ import annotations
import argparse
import os
import sys
from itertools import zip_longest
from typing import Optional

import fastb
from fastb.alphabet import AlphabetError
from fastb.io import (iter_fasta, iter_records, open_text_or_gz, stream_fasta, write_records,
                      wrap_bytes)


def _out():
    return sys.stdout.buffer


def _emit(out, name: str, seq: str, wrap: int) -> None:
    out.write(b">" + name.encode("utf-8") + b"\n" + wrap_bytes(seq, wrap))


def _wrap_arg(p):
    p.add_argument("-w", "--wrap", type=int, default=80, metavar="N",
                   help="Wrap sequence to N chars (0=no wrap, default 80).")


def _threads_arg(p):
    p.add_argument("-p", "--threads", type=int, default=1, metavar="N",
                   help="Decode chunks in N worker processes (default 1).")


def _strip_ext(path: str, exts) -> str:
    low = path.lower()
    for ext in exts:
        if low.endswith(ext):
            return path[:-len(ext)]
    return path


# ---------------------------------------------------------------------------
# cat / view
# ---------------------------------------------------------------------------

def _add_cat_parser(sub, name, help_text):
    p = sub.add_parser(name, help=help_text)
    p.add_argument("file", metavar="file.fastb")
    p.add_argument("-n", "--names", metavar="NAME[,NAME...]",
                   help="Only emit these records (comma-separated).")
    _wrap_arg(p)
    _threads_arg(p)
    p.set_defaults(func=_cmd_cat, cmd=name)


def _cmd_cat(args) -> int:
    names = set(args.names.split(",")) if args.names else None
    try:
        stream_fasta(args.file, _out(), args.wrap, args.threads, names)
    except ValueError as e:
        print(f"fastb {args.cmd}: {e}", file=sys.stderr)
        return 3
    if names:
        missing = names - set(fastb.read_index(args.file))
        for name in sorted(missing):
            print(f"fastb {args.cmd}: record {name!r} not found", file=sys.stderr)
        return 2 if missing else 0
    return 0


# ---------------------------------------------------------------------------
# head
# ---------------------------------------------------------------------------

def _add_head_parser(sub):
    p = sub.add_parser("head", help="Print the first N records as FASTA.")
    p.add_argument("file", metavar="file.fastb")
    p.add_argument("-n", type=int, default=5, metavar="N", help="Number of records (default 5).")
    _wrap_arg(p)
    p.set_defaults(func=_cmd_head)


def _cmd_head(args) -> int:
    out = _out()
    try:
        for i, (name, seq) in enumerate(iter_records(args.file)):
            if i >= args.n:
                break
            _emit(out, name, seq, args.wrap)
    except ValueError as e:
        print(f"fastb head: {e}", file=sys.stderr)
        return 3
    return 0


# ---------------------------------------------------------------------------
# stats
# ---------------------------------------------------------------------------

def _add_stats_parser(sub):
    p = sub.add_parser("stats", help="Per-record statistics.")
    p.add_argument("file", metavar="file.fastb")
    p.set_defaults(func=_cmd_stats)


def _cmd_stats(args) -> int:
    import numpy as np
    print("name\tlength\talpha\tenc\tgc_pct\tn_pct\tmasked_pct")
    tot = np.zeros(4, dtype=np.int64)  # length, gc, n, masked
    try:
        with open(args.file, "rb") as f:
            for rec in fastb.read_file(f):
                enc, *_ = fastb.encode_sequence(rec.sequence, rec.alpha)
                counts = fastb.alphabet._letter_counts(rec.sequence)
                raw = np.frombuffer(rec.sequence.encode("ascii"), dtype=np.uint8)
                row = np.array([len(raw), counts[ord("G")] + counts[ord("C")],
                                counts[ord("N")], int(((raw >= 97) & (raw <= 122)).sum())])
                tot += row
                pct = [f"{100 * x / row[0]:.2f}" if row[0] else "0.00" for x in row[1:]]
                print(f"{rec.name}\t{row[0]}\t{rec.alpha}\t{enc}\t" + "\t".join(pct))
    except ValueError as e:
        print(f"fastb stats: {e}", file=sys.stderr)
        return 3
    if tot[0]:
        pct = [f"{100 * x / tot[0]:.2f}" for x in tot[1:]]
        print(f"TOTAL\t{tot[0]}\t-\t-\t" + "\t".join(pct))
    return 0


# ---------------------------------------------------------------------------
# extract
# ---------------------------------------------------------------------------

def _add_extract_parser(sub):
    p = sub.add_parser("extract", help="Random-access extraction to FASTA.")
    p.add_argument("file", metavar="file.fastb")
    p.add_argument("names", nargs="+", metavar="NAME")
    _wrap_arg(p)
    p.add_argument("--region", metavar="NAME:START-END",
                   help="1-based inclusive subrange of one record.")
    p.set_defaults(func=_cmd_extract)


def _cmd_extract(args) -> int:
    out = _out()
    try:
        idx = fastb.read_index(args.file)
    except ValueError as e:
        print(f"fastb extract: {e}", file=sys.stderr)
        return 3

    region = None
    if args.region:
        try:
            rname, rrange = args.region.rsplit(":", 1)
            rstart, rend = rrange.split("-")
            region = (rname, int(rstart) - 1, int(rend))
        except ValueError:
            print("fastb extract: invalid --region format, expected NAME:START-END",
                  file=sys.stderr)
            return 1

    missing = [n for n in args.names if n not in idx]
    for name in missing:
        print(f"fastb extract: record {name!r} not found", file=sys.stderr)

    for name in args.names:
        if name in missing:
            continue
        rec = fastb.read_record(args.file, name)
        seq, display = rec.sequence, rec.name
        if region and region[0] == name:
            seq = seq[region[1]:region[2]]
            display = f"{rec.name}:{region[1] + 1}-{region[2]}"
        _emit(out, display, seq, args.wrap)
    return 2 if missing else 0


# ---------------------------------------------------------------------------
# encode / decode
# ---------------------------------------------------------------------------

def _add_encode_parser(sub):
    p = sub.add_parser("encode", help="Convert FASTA (.fasta, .fa, .fna, or .gz) to FASTB. "
                                      "Several files are encoded in one process, each on its own: "
                                      "a failed file leaves no output and the exit code is the "
                                      "highest per-file code.")
    p.add_argument("files", nargs="+", metavar="file.fasta")
    p.add_argument("-o", "--output", metavar="PATH",
                   help="Output path, one input only (default: input with the extension "
                        "replaced by .fastb).")
    p.add_argument("--append", action="store_true",
                   help="Name each output <input>.fastb, keeping the input's extension, so "
                        "x.fa and x.fasta cannot collide.")
    p.add_argument("--comment", metavar="TEXT", help="Optional file-level comment.")
    p.add_argument("--force-nucleotide", action="store_true",
                   help="Skip the protein heuristic. Never bypasses the hard F/I/L/P/Q/E/Z reject.")
    p.add_argument("--verify", action="store_true",
                   help="After writing, decode the output and compare it with the input record by "
                        "record (full header line and sequence). On any difference, including a "
                        "dropped header description, remove the output and exit 3.")
    p.set_defaults(func=_cmd_encode)


def _fasta_lines(path: str):
    """Yield (header, sequence) bytes per record, header without '>' and with its description."""
    header, parts = None, []
    with open_text_or_gz(path) as f:
        for line in f:
            line = line.rstrip(b"\r\n")
            if line.startswith(b">"):
                if header is not None:
                    yield header, b"".join(parts)
                header, parts = line[1:], []
            elif line:
                parts.append(line)
    if header is not None:
        yield header, b"".join(parts)


def _first_difference(fasta_path: str, fastb_path: str) -> Optional[str]:
    """Return what differs between the FASTA and the written FASTB, or None if nothing does."""
    pairs = zip_longest(_fasta_lines(fasta_path), iter_records(fastb_path))
    for i, (a, b) in enumerate(pairs, 1):
        if a is None or b is None:
            return f"record counts differ at record {i}"
        if a[0] != b[0].encode("utf-8"):
            return f"record {i}: header '{a[0].decode('utf-8', 'replace')}' became '{b[0]}'"
        if a[1] != b[1].encode("ascii"):
            return f"record {i} ({b[0]}): sequence differs"
    return None


def _encode_one(in_path: str, out_path: str, args) -> int:
    try:
        n = write_records(out_path, iter_fasta(in_path, as_bytes=True), file_comment=args.comment,
                          force_nucleotide=args.force_nucleotide)
    except AlphabetError as e:
        print(f"fastb encode: {in_path}: {e}", file=sys.stderr)
        return 4
    except (OSError, ValueError, UnicodeDecodeError) as e:
        print(f"fastb encode: {in_path}: {e}", file=sys.stderr)
        return 1
    if args.verify:
        diff = _first_difference(in_path, out_path)
        if diff:
            os.remove(out_path)
            print(f"fastb encode: {in_path}: verify failed, {diff}; {out_path} removed",
                  file=sys.stderr)
            return 3
    print(f"Encoded {n} record(s) -> {out_path}", file=sys.stderr)
    return 0


def _cmd_encode(args) -> int:
    if args.output and (len(args.files) > 1 or args.append):
        print("fastb encode: -o takes one input and cannot be combined with --append",
              file=sys.stderr)
        return 1
    rc = 0
    for in_path in args.files:
        if args.output:
            out_path = args.output
        elif args.append:
            out_path = in_path + ".fastb"
        else:
            out_path = _strip_ext(_strip_ext(in_path, (".gz",)), (".fasta", ".fa", ".fna")) + ".fastb"
        rc = max(rc, _encode_one(in_path, out_path, args))
    return rc


def _add_decode_parser(sub):
    p = sub.add_parser("decode", help="Convert FASTB to FASTA.")
    p.add_argument("file", metavar="file.fastb")
    p.add_argument("-o", "--output", metavar="PATH",
                   help="Output path (default: replace .fastb with .fasta).")
    _wrap_arg(p)
    _threads_arg(p)
    p.set_defaults(func=_cmd_decode)


def _cmd_decode(args) -> int:
    out_path = args.output or _strip_ext(args.file, (".fastb",)) + ".fasta"
    try:
        with open(out_path, "wb") as fout:
            stream_fasta(args.file, fout, args.wrap, args.threads)
    except ValueError as e:
        print(f"fastb decode: {e}", file=sys.stderr)
        return 3
    except OSError as e:
        print(f"fastb decode: {e}", file=sys.stderr)
        return 1
    print(f"Decoded -> {out_path}", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# index / info / verify
# ---------------------------------------------------------------------------

def _add_index_parser(sub):
    p = sub.add_parser("index", help="Print the record index: number, name, offset, length.")
    p.add_argument("file", metavar="file.fastb")
    p.set_defaults(func=_cmd_index)


def _cmd_index(args) -> int:
    try:
        idx = fastb.read_index(args.file)
    except ValueError as e:
        print(f"fastb index: {e}", file=sys.stderr)
        return 3
    print("record\tname\toffset\tlength")
    for k, (name, (off, length)) in enumerate(idx.items()):
        print(f"{k}\t{name}\t{off}\t{length}")
    return 0


def _add_info_parser(sub):
    p = sub.add_parser("info", help="File summary: records, bytes, bits per base, tier counts.")
    p.add_argument("file", metavar="file.fastb")
    p.set_defaults(func=_cmd_info)


def _cmd_info(args) -> int:
    records = bases = enc2 = enc4 = nruns = mask = iupac = 0
    try:
        with open(args.file, "rb") as f:
            for fields, _, _ in fastb.iter_raw(f, read_payload=False):
                records += 1
                bases += int(fields["LEN"])
                if fields["ENC"] == "2":
                    enc2 += 1
                else:
                    enc4 += 1
                nruns += fields.get("NRUNS", "").count(":")
                mask += fields.get("MASK", "").count(":")
                iupac += len(fields["IUPAC"].split(",")) if fields.get("IUPAC") else 0
    except ValueError as e:
        print(f"fastb info: {e}", file=sys.stderr)
        return 3
    size = os.path.getsize(args.file)
    print(f"file\t{args.file}")
    print(f"records\t{records}")
    print(f"bases\t{bases}")
    print(f"bytes\t{size}")
    print(f"bits_per_base\t{8 * size / bases:.3f}" if bases else "bits_per_base\tn/a")
    print(f"records_2bit\t{enc2}")
    print(f"records_4bit\t{enc4}")
    print(f"n_runs\t{nruns}")
    print(f"iupac_runs\t{iupac}")
    print(f"mask_runs\t{mask}")
    return 0


def _add_verify_parser(sub):
    p = sub.add_parser("verify", help="Check magic, CRCs, index, and alphabet compliance.")
    p.add_argument("file", metavar="file.fastb")
    p.set_defaults(func=_cmd_verify)


def _cmd_verify(args) -> int:
    errors = []
    count = 0
    try:
        with open(args.file, "rb") as f:
            for _ in fastb.read_file(f):
                count += 1
    except ValueError as e:
        errors.append(str(e))
    if not errors:
        try:
            idx = fastb.read_index(args.file)
            if len(idx) != count:
                errors.append(f"index lists {len(idx)} records, file has {count}")
            for name in idx:
                try:
                    fastb.read_record(args.file, name)
                except (KeyError, ValueError) as e:
                    errors.append(f"index inconsistency for {name!r}: {e}")
        except ValueError as e:
            errors.append(f"index: {e}")
    if errors:
        for err in errors:
            print(f"ERROR: {err}", file=sys.stderr)
        print(f"FAIL: {len(errors)} error(s) in {args.file}")
        return 3
    print(f"OK: {count} record(s), index consistent: {args.file}")
    return 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="fastb", description="FASTB: binary sequence format.")
    parser.add_argument("--version", action="version", version=f"fastb {fastb.__version__}")
    sub = parser.add_subparsers(dest="command", required=True)
    _add_encode_parser(sub)
    _add_decode_parser(sub)
    _add_cat_parser(sub, "cat", "Write records as FASTA to stdout (for pipes).")
    _add_cat_parser(sub, "view", "Alias of cat.")
    _add_head_parser(sub)
    _add_extract_parser(sub)
    _add_index_parser(sub)
    _add_info_parser(sub)
    _add_stats_parser(sub)
    _add_verify_parser(sub)
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except BrokenPipeError:
        return 0
    except Exception:
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
