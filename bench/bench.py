"""Measure fastb against pigz/gzip on one FASTA file.

Usage: python bench/bench.py genome.fasta [-p THREADS]

Reports bytes on disk, encode and decode wall time, and peak RSS for both
tools. pigz is used when found on PATH, else gzip (single thread). Peak RSS
needs the `psutil` package (pip install fastb[bench]).
"""

from __future__ import annotations
import argparse
import os
import shutil
import subprocess
import sys
import time


def _peak_rss(cmd) -> tuple[float, int]:
    """Run cmd, return (wall seconds, peak RSS in bytes or -1 if psutil missing)."""
    try:
        import psutil
    except ImportError:
        t = time.perf_counter()
        subprocess.run(cmd, check=True)
        return time.perf_counter() - t, -1
    t = time.perf_counter()
    proc = psutil.Popen(cmd)
    peak = 0
    while proc.poll() is None:
        try:
            peak = max(peak, proc.memory_info().rss)
        except psutil.NoSuchProcess:
            break
        time.sleep(0.01)
    if proc.returncode:
        raise SystemExit(f"{cmd[0]} failed with exit code {proc.returncode}")
    return time.perf_counter() - t, peak


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("-p", "--threads", type=int, default=1)
    args = ap.parse_args()

    fasta = args.fasta
    work = fasta + ".bench"
    os.makedirs(work, exist_ok=True)
    copy = os.path.join(work, "in.fasta")
    shutil.copyfile(fasta, copy)

    gz_tool = "pigz" if shutil.which("pigz") else "gzip"
    py = sys.executable
    rows = []

    # gzip / pigz
    gz_args = ["-6", "-f", "-k"] + (["-p", str(args.threads)] if gz_tool == "pigz" else [])
    t_enc, m_enc = _peak_rss([gz_tool] + gz_args + [copy])
    gz_size = os.path.getsize(copy + ".gz")
    dec_cmd = [gz_tool, "-dc"] + (["-p", str(args.threads)] if gz_tool == "pigz" else [])
    out = open(os.path.join(work, "gz.out"), "wb")
    t = time.perf_counter()
    subprocess.run(dec_cmd + [copy + ".gz"], stdout=out, check=True)
    t_dec = time.perf_counter() - t
    out.close()
    rows.append((f"{gz_tool} -6", gz_size, t_enc, t_dec, m_enc))

    # fastb
    fb = os.path.join(work, "in.fastb")
    t_enc, m_enc = _peak_rss([py, "-m", "fastb.cli", "encode", copy, "-o", fb])
    fb_size = os.path.getsize(fb)
    out = open(os.path.join(work, "fb.out"), "wb")
    t = time.perf_counter()
    subprocess.run([py, "-m", "fastb.cli", "cat", fb], stdout=out, check=True)
    t_dec = time.perf_counter() - t
    out.close()
    rows.append(("fastb", fb_size, t_enc, t_dec, m_enc))

    fa_size = os.path.getsize(copy)
    print(f"input: {fasta} ({fa_size / 1e6:.2f} MB), threads={args.threads}")
    print("| tool | bytes on disk | encode s | decode s | peak RSS MB |")
    print("|---|---|---|---|---|")
    for name, size, te, td, rss in rows:
        rss_s = f"{rss / 1e6:.0f}" if rss >= 0 else "n/a"
        print(f"| {name} | {size:,} | {te:.2f} | {td:.2f} | {rss_s} |")
    shutil.rmtree(work)


if __name__ == "__main__":
    main()
