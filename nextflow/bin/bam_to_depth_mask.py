#!/usr/bin/env python3
"""Emit contiguous callable BED intervals using per-base, quality-filtered BAM depth."""

import argparse
import subprocess


def make_mask(bam, output, minimum=2, maximum=20):
    if minimum < 1 or maximum < minimum:
        raise ValueError("Depth bounds must satisfy 1 <= min_depth <= max_depth")
    process = subprocess.Popen(["samtools", "depth", "-a", "-q", "20", "-Q", "20", bam],
                               stdout=subprocess.PIPE, text=True)
    intervals = 0
    with open(output, "w") as handle:
        chrom, start, end = None, None, None
        for line in process.stdout:
            current, position, depth = line.split()
            position, depth = int(position), int(depth)
            if minimum <= depth <= maximum:
                if current == chrom and position - 1 == end:
                    end = position
                    continue
                if chrom is not None:
                    handle.write(f"{chrom}\t{start}\t{end}\n")
                    intervals += 1
                chrom, start, end = current, position - 1, position
            elif chrom is not None:
                handle.write(f"{chrom}\t{start}\t{end}\n")
                intervals += 1
                chrom, start, end = None, None, None
        if chrom is not None:
            handle.write(f"{chrom}\t{start}\t{end}\n")
            intervals += 1
    result = process.wait()
    if result != 0:
        raise RuntimeError(f"samtools depth exited with status {result}")
    if intervals == 0:
        raise ValueError(f"No callable bases at depths {minimum}..{maximum} in {bam}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam")
    parser.add_argument("out_bed")
    parser.add_argument("--min-depth", type=int, default=2)
    parser.add_argument("--max-depth", type=int, default=20)
    args = parser.parse_args()
    make_mask(args.bam, args.out_bed, args.min_depth, args.max_depth)


if __name__ == "__main__":
    main()
