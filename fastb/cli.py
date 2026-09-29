"""FASTB v3 command-line interface.

Single entry point `fastb` with subcommand dispatch.
Exit codes: 0=success, 1=generic error, 2=partial success,
            3=file format error, 4=alphabet error.
"""

from __future__ import annotations
import argparse
import sys

import fastb
from fastb.alphabet import AlphabetError, ProteinDetectedError


def _wrap(seq: str, width: int) -> str:
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def _emit_fasta(name: str, seq: str, wrap: int) -> None:
    sys.stdout.write(f">{name}\n{_wrap(seq, wrap)}\n")


# ---------------------------------------------------------------------------
# view
# ---------------------------------------------------------------------------

def _add_view_parser(subparsers):
    p = subparsers.add_parser("view", help="Decode records as FASTA to stdout.")
    p.add_argument("file", metavar="file.fastb")
    p.add_argument("-n", "--names", metavar="NAME[,NAME...]",
                   help="Only emit these records (comma-separated).")
    p.add_argument("-w", "--wrap", type=int, default=80, metavar="N",
                   help="Wrap sequence to N chars (0=no wrap, default 80).")
    p.add_argument("--no-index", action="store_true",
                   help="Stream from the top; ignore footer index.")
    p.set_defaults(func=_cmd_view)


def _cmd_view(args) -> int:
    names = set(args.names.split(",")) if args.names else None
    wrap = args.wrap

    if names and not args.no_index:
        try:
            idx = fastb.read_index(args.file)
        except Exception:
            idx = None
        if idx is not None:
            missing = names - set(idx)
            for name in sorted(missing):
                print(f"fastb view: record {name!r} not found", file=sys.stderr)
            for name in names - missing:
                try:
                    rec = fastb.read_record(args.file, name)
                    _emit_fasta(rec.name, rec.sequence, wrap)
                except (KeyError, ValueError) as e:
                    print(f"fastb view: {e}", file=sys.stderr)
                    return 3
            return 2 if missing else 0

    try:
        with open(args.file, "rb") as f:
            for rec in fastb.read_file(f):
                if names is None or rec.name in names:
                    _emit_fasta(rec.name, rec.sequence, wrap)
    except ValueError as e:
        print(f"fastb view: {e}", file=sys.stderr)
        return 3
    return 0


# ---------------------------------------------------------------------------
# head
# ---------------------------------------------------------------------------

def _add_head_parser(subparsers):
    p = subparsers.add_parser("head", help="Print the first N records as FASTA.")
    p.add_argument("file", metavar="file.fastb")
    p.add_argument("-n", type=int, default=5, metavar="N",
                   help="Number of records (default 5).")
    p.add_argument("--bases", type=int, default=0, metavar="N",
                   help="Stop once ~N bases have been emitted.")
    p.add_argument("-w", "--wrap", type=int, default=80, metavar="N",
                   help="Wrap sequence to N chars (0=no wrap, default 80).")
    p.set_defaults(func=_cmd_head)


def _cmd_head(args) -> int:
    count = 0
    bases_emitted = 0
    try:
        with open(args.file, "rb") as f:
            for rec in fastb.read_file(f):
                if count >= args.n:
                    break
                if args.bases and bases_emitted >= args.bases:
                    break
                _emit_fasta(rec.name, rec.sequence, args.wrap)
                count += 1
                bases_emitted += len(rec.sequence)
    except ValueError as e:
        print(f"fastb head: {e}", file=sys.stderr)
        return 3
    return 0


# ---------------------------------------------------------------------------
# stats
# ---------------------------------------------------------------------------

def _add_stats_parser(subparsers):
    p = subparsers.add_parser("stats", help="Per-record statistics.")
    p.add_argument("file", metavar="file.fastb")
    p.set_defaults(func=_cmd_stats)


def _cmd_stats(args) -> int:
    print("name\tlength\talpha\tenc\tgc_pct\tn_pct\tmasked_pct")
    total_len = total_gc = total_n = total_masked = count = 0
    try:
        with open(args.file, "rb") as f:
            for rec in fastb.read_file(f):
                seq_up = rec.sequence.upper()
                length = len(seq_up)
                basic = "ACGT" if rec.alpha == "D" else "ACGU"
                enc = 2 if all(c in basic for c in seq_up) else 4
                gc = sum(1 for c in seq_up if c in "GC")
                n = seq_up.count("N")
                masked = sum(1 for c in rec.sequence if c.islower())
                def pct(x):
                    return f"{100 * x / length:.2f}" if length else "0.00"
                print(f"{rec.name}\t{length}\t{rec.alpha}\t{enc}\t"
                      f"{pct(gc)}\t{pct(n)}\t{pct(masked)}")
                total_len += length
                total_gc += gc
                total_n += n
                total_masked += masked
                count += 1
    except ValueError as e:
        print(f"fastb stats: {e}", file=sys.stderr)
        return 3

    if count:
        def tot_pct(x):
            return f"{100 * x / total_len:.2f}" if total_len else "0.00"
        print(f"TOTAL\t{total_len}\t-\t-\t"
              f"{tot_pct(total_gc)}\t{tot_pct(total_n)}\t{tot_pct(total_masked)}")
    return 0


# ---------------------------------------------------------------------------
# extract
# ---------------------------------------------------------------------------

def _add_extract_parser(subparsers):
    p = subparsers.add_parser("extract", help="Random-access extraction to FASTA.")
    p.add_argument("file", metavar="file.fastb")
    p.add_argument("names", nargs="+", metavar="NAME")
    p.add_argument("-w", "--wrap", type=int, default=80, metavar="N",
                   help="Wrap sequence to N chars (0=no wrap, default 80).")
    p.add_argument("--region", metavar="NAME:START-END",
                   help="1-based inclusive subrange of one record.")
    p.set_defaults(func=_cmd_extract)


def _cmd_extract(args) -> int:
    try:
        idx = fastb.read_index(args.file)
    except ValueError as e:
        print(f"fastb extract: {e}", file=sys.stderr)
        return 3

    region_name = region_start = region_end = None
    if args.region:
        try:
            rname, rrange = args.region.rsplit(":", 1)
            rstart, rend = rrange.split("-")
            region_name = rname
            region_start = int(rstart) - 1  # 0-based
            region_end = int(rend)           # end is inclusive, slice [start:end]
        except (ValueError, AttributeError):
            print("fastb extract: invalid --region format, expected NAME:START-END",
                  file=sys.stderr)
            return 1

    missing = [n for n in args.names if n not in idx]
    for name in missing:
        print(f"fastb extract: record {name!r} not found", file=sys.stderr)

    exit_code = 2 if missing else 0
    for name in args.names:
        if name in missing:
            continue
        try:
            rec = fastb.read_record(args.file, name)
        except (KeyError, ValueError) as e:
            print(f"fastb extract: {e}", file=sys.stderr)
            exit_code = 2
            continue

        seq = rec.sequence
        display_name = rec.name
        if region_name == name and region_start is not None:
            seq = seq[region_start:region_end]
            display_name = f"{rec.name}:{region_start + 1}-{region_end}"

        _emit_fasta(display_name, seq, args.wrap)

    return exit_code


# ---------------------------------------------------------------------------
# encode
# ---------------------------------------------------------------------------

def _add_encode_parser(subparsers):
    p = subparsers.add_parser("encode", help="Convert FASTA to FASTB v3.")
    p.add_argument("file", metavar="file.fasta")
    p.add_argument("-o", "--output", metavar="PATH",
                   help="Output path (default: replace .fasta/.fa with .fastb).")
    p.add_argument("--comment", metavar="TEXT",
                   help="Optional file-level comment.")
    p.add_argument("--force-nucleotide", action="store_true",
                   help="Skip the protein-detection heuristic (rule 3). "
                        "Does NOT bypass the hard F/I/L/P/Q/E/Z reject.")
    p.set_defaults(func=_cmd_encode)


def _cmd_encode(args) -> int:
    from fastb.alphabet import require_nucleotide

    out_path = args.output
    if not out_path:
        base = args.file
        for ext in (".fasta", ".fa", ".FASTA", ".FA"):
            if base.endswith(ext):
                base = base[:-len(ext)]
                break
        out_path = base + ".fastb"

    try:
        raw_records = list(_parse_fasta(args.file))
    except OSError as e:
        print(f"fastb encode: {e}", file=sys.stderr)
        return 1

    # Validate all records before writing anything (no partial files)
    fastb_records = []
    for name, seq in raw_records:
        seq_upper = seq.upper()
        has_t = "T" in seq_upper
        has_u = "U" in seq_upper

        if has_t and has_u:
            print(f"fastb encode: record {name!r} contains both T and U; "
                  "cannot determine alphabet.", file=sys.stderr)
            return 4

        if has_u and not has_t:
            alpha = "R"
        elif has_t and not has_u:
            alpha = "D"
        else:
            print(f"fastb encode: record {name!r} contains neither T nor U; "
                  "defaulting to ALPHA=D.", file=sys.stderr)
            alpha = "D"

        try:
            require_nucleotide(seq_upper, alpha, name,
                               force_nucleotide=args.force_nucleotide)
        except ProteinDetectedError as e:
            print(f"fastb encode: {e}", file=sys.stderr)
            return 4
        except AlphabetError as e:
            print(f"fastb encode: {e}", file=sys.stderr)
            return 4

        fastb_records.append(fastb.Record(name, seq, alpha=alpha))

    # All records validated. Now write atomically via in-memory buffer
    import io as _io
    buf = _io.BytesIO()
    fastb.write_file(fastb_records, buf, file_comment=args.comment)

    try:
        with open(out_path, "wb") as f:
            f.write(buf.getvalue())
    except OSError as e:
        print(f"fastb encode: {e}", file=sys.stderr)
        return 1

    print(f"Encoded {len(fastb_records)} record(s) \u2192 {out_path}", file=sys.stderr)
    return 0


def _parse_fasta(path: str):
    """Yield (name, sequence) pairs from a FASTA file. Stdlib only."""
    name = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                header = line[1:].split()
                name = header[0] if header else ""
                parts = []
            elif line:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts)


# ---------------------------------------------------------------------------
# decode
# ---------------------------------------------------------------------------

def _add_decode_parser(subparsers):
    p = subparsers.add_parser("decode", help="Convert FASTB v3 to FASTA.")
    p.add_argument("file", metavar="file.fastb")
    p.add_argument("-o", "--output", metavar="PATH",
                   help="Output path (default: replace .fastb with .fasta).")
    p.add_argument("-w", "--wrap", type=int, default=80, metavar="N",
                   help="Wrap sequence to N chars (0=no wrap, default 80).")
    p.set_defaults(func=_cmd_decode)


def _cmd_decode(args) -> int:
    out_path = args.output
    if not out_path:
        base = args.file
        if base.endswith(".fastb"):
            base = base[:-6]
        out_path = base + ".fasta"

    try:
        with open(args.file, "rb") as fin, open(out_path, "w") as fout:
            for rec in fastb.read_file(fin):
                fout.write(f">{rec.name}\n{_wrap(rec.sequence, args.wrap)}\n")
    except ValueError as e:
        print(f"fastb decode: {e}", file=sys.stderr)
        return 3
    except OSError as e:
        print(f"fastb decode: {e}", file=sys.stderr)
        return 1

    print(f"Decoded \u2192 {out_path}", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# verify
# ---------------------------------------------------------------------------

def _add_verify_parser(subparsers):
    p = subparsers.add_parser("verify",
                              help="Check magic, CRCs, index, and alphabet compliance.")
    p.add_argument("file", metavar="file.fastb")
    p.set_defaults(func=_cmd_verify)


def _cmd_verify(args) -> int:
    errors = []
    record_count = 0

    try:
        with open(args.file, "rb") as f:
            for rec in fastb.read_file(f):
                record_count += 1
    except ValueError as e:
        errors.append(str(e))

    if not errors:
        try:
            idx = fastb.read_index(args.file)
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

    print(f"OK: {record_count} record(s), index consistent \u2014 {args.file}")
    return 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="fastb",
        description="FASTB v3: binary sequence format.",
    )
    parser.add_argument(
        "--version", action="version", version="fastb 3.0.0"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    _add_view_parser(subparsers)
    _add_head_parser(subparsers)
    _add_stats_parser(subparsers)
    _add_extract_parser(subparsers)
    _add_encode_parser(subparsers)
    _add_decode_parser(subparsers)
    _add_verify_parser(subparsers)

    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except Exception as e:
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
