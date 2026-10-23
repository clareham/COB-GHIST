"""Validated BED intervals and callable-base counting for VCF converters."""

import bisect
from collections import defaultdict


def load_mask(path, allow_empty=False):
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
    if not result and not allow_empty:
        raise ValueError(f"{path}: empty callability mask")
    return result


def full_mask(lengths):
    """Callable mask covering every declared contig end to end."""
    return {chrom: ([0], [length]) for chrom, length in lengths.items()}


def subtract_mask(mask, excluded):
    """Remove `excluded` intervals (same {chrom: (starts, ends)} shape) from `mask`."""
    result = {}
    for chrom, (starts, ends) in mask.items():
        cut_starts, cut_ends = excluded.get(chrom, ([], []))
        kept = []
        index = 0
        for start, end in zip(starts, ends):
            cursor = start
            while index < len(cut_starts) and cut_ends[index] <= cursor:
                index += 1
            probe = index
            while probe < len(cut_starts) and cut_starts[probe] < end:
                if cut_starts[probe] > cursor:
                    kept.append((cursor, cut_starts[probe]))
                cursor = max(cursor, cut_ends[probe])
                probe += 1
            if cursor < end:
                kept.append((cursor, end))
        if kept:
            result[chrom] = ([s for s, _ in kept], [e for _, e in kept])
    return result


def resolve_mask(lengths, mask_path=None, exclude_path=None, assume_callable=False):
    """Build the callable mask from an optional BED, assumed-callable contigs and exclusions.

    Returns None when everything is callable and nothing is excluded, so
    callers can keep their fast "no mask" path.
    """
    if mask_path is None and not assume_callable:
        raise ValueError("Provide a callability BED or explicitly choose --assume-callable")
    mask = load_mask(mask_path) if mask_path else None
    if exclude_path:
        excluded = load_mask(exclude_path, allow_empty=True)
        if excluded:
            mask = subtract_mask(mask if mask is not None else full_mask(lengths), excluded)
    return mask


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


def iter_positions(records, chrom_of=lambda rec: rec.CHROM, pos_of=lambda rec: rec.POS):
    """Group consecutive records sharing CHROM/POS.

    Simulated GHIST VCFs stack multiple mutations at one coordinate as
    separate records; callers treat such groups like multiallelic sites.
    Ordering within each contig is validated here so that callers only
    have to handle complete per-position groups.
    """
    group = []
    previous = {}
    for rec in records:
        chrom, pos = chrom_of(rec), pos_of(rec)
        if group and (chrom != chrom_of(group[0]) or pos != pos_of(group[0])):
            yield group
            group = []
        if not group:
            if pos <= previous.get(chrom, 0):
                raise ValueError(f"VCF positions must be sorted within each contig: {chrom}:{pos}")
            previous[chrom] = pos
        group.append(rec)
    if group:
        yield group


def contig_lengths(vcf):
    lengths = dict(zip(vcf.seqnames, vcf.seqlens))
    if not lengths or any(length is None or length <= 0 for length in lengths.values()):
        raise ValueError("VCF must declare a positive length for every contig")
    return lengths
