#!/usr/bin/env python3
"""Build a simple per-window callability mask BED from one BAM's per-base
depth, for low-coverage data where no mask is otherwise provided.

msmc-tools' own `bamCaller.py` does this alongside VCF generation, but
it isn't packaged on bioconda/PyPI (same situation as
generate_multihetsep.py - see bin/vcf_to_multihetsep.py), so this is a
minimal standalone reimplementation: a window (default 100bp, matching
PSMC/MSMC2 convention) is "callable" if its mean per-base depth falls
within [min_depth, max_depth]. Pass the resulting BED to
bin/vcf_to_multihetsep.py's --mask so MSMC2's `nSitesSincePrevCall`
counts reflect actual low-coverage callability instead of raw physical
distance.
"""
import argparse
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("bam")
    ap.add_argument("out_bed")
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--min-depth", type=int, default=2)
    ap.add_argument("--max-depth", type=int, default=20)
    args = ap.parse_args()

    proc = subprocess.Popen(
        ["samtools", "depth", "-a", args.bam],
        stdout=subprocess.PIPE, text=True,
    )

    cur_chrom, win_start, depths = None, None, []

    def flush(out):
        if cur_chrom is None or not depths:
            return
        mean_depth = sum(depths) / len(depths)
        if args.min_depth <= mean_depth <= args.max_depth:
            out.write(f"{cur_chrom}\t{win_start}\t{win_start + len(depths)}\n")

    with open(args.out_bed, "w") as out:
        for line in proc.stdout:
            chrom, pos, depth = line.rstrip("\n").split("\t")
            pos, depth = int(pos), int(depth)
            win = (pos - 1) // args.window
            this_start = win * args.window
            if chrom != cur_chrom or this_start != win_start:
                flush(out)
                cur_chrom, win_start, depths = chrom, this_start, []
            depths.append(depth)
        flush(out)

    ret = proc.wait()
    if ret != 0:
        sys.exit(f"samtools depth exited with status {ret}")


if __name__ == "__main__":
    main()
