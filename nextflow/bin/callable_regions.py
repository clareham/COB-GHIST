"""Validated BED intervals and callable-base counting for VCF converters."""

import bisect
from collections import defaultdict


def load_mask(path):
    by_chrom = defaultdict(list)
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{number}: expected chromosome, start, end")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{number}: invalid BED interval")
            by_chrom[chrom].append((start, end))
    result = {}
    for chrom, intervals in by_chrom.items():
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        result[chrom] = ([start for start, _ in merged], [end for _, end in merged])
    if not result:
        raise ValueError(f"{path}: empty callability mask")
    return result


def callable_bases_in_range(intervals, start_excl, end_incl):
    """Count BED-covered bases in the 1-based interval (start, end]."""
    if end_incl < start_excl:
        raise ValueError("Positions must increase within each contig")
    if intervals is None:
        return end_incl - start_excl
    starts, ends = intervals
    index = bisect.bisect_right(ends, start_excl)
    total = 0
    while index < len(starts) and starts[index] < end_incl:
        total += min(ends[index], end_incl) - max(starts[index], start_excl)
        index += 1
    return total


def is_callable(intervals, pos):
    return callable_bases_in_range(intervals, pos - 1, pos) == 1


def contig_lengths(vcf):
    lengths = dict(zip(vcf.seqnames, vcf.seqlens))
    if not lengths or any(length is None or length <= 0 for length in lengths.values()):
        raise ValueError("VCF must declare a positive length for every contig")
    return lengths
