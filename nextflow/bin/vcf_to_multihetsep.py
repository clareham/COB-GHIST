#!/usr/bin/env python3
"""Create MSMC2 multihetsep input, preserving unphased genotype ambiguity.

An explicit --assume-callable is required without a BED. Invariant sites
are not emitted, but callable invariant bases contribute to distances.
Use --split to write one input file per contig (required by MSMC2).
"""

import argparse
import itertools
import sys
from contextlib import ExitStack
from pathlib import Path

from cyvcf2 import VCF

from callable_regions import callable_bases_in_range, contig_lengths, is_callable, load_mask


def convert(vcf_path, output, samples=None, mask_path=None, assume_callable=False, split=False):
    if mask_path is None and not assume_callable:
        raise ValueError("Provide a callability BED or explicitly choose --assume-callable")
    vcf = VCF(vcf_path)
    selected = list(samples) if samples is not None else vcf.samples[:2]
    if not selected or len(set(selected)) != len(selected) or set(selected) - set(vcf.samples):
        raise ValueError("Requested samples must be a nonempty, unique subset of VCF samples")
    vcf.set_samples(selected)
    lengths = contig_lengths(vcf)
    mask = load_mask(mask_path) if mask_path else None
    previous_record = {}
    previous_emitted = {}
    missing_positions = {}
    excluded_until = {}
    handles = {}
    paths = []
    skipped = 0
    with ExitStack() as stack:
        if split:
            Path(output).mkdir(parents=True, exist_ok=True)
        else:
            handles[None] = stack.enter_context(open(output, "w"))
        for rec in vcf:
            chrom, pos = rec.CHROM, rec.POS
            if chrom not in lengths or not 1 <= pos <= lengths[chrom]:
                raise ValueError(f"Record outside declared contig: {chrom}:{pos}")
            if pos <= previous_record.get(chrom, 0):
                raise ValueError("Positions must be strictly increasing within each contig")
            previous_record[chrom] = pos
            intervals = mask.get(chrom, ([], [])) if mask is not None else None
            end = min(pos + len(rec.REF) - 1, lengths[chrom])
            blocked_end = excluded_until.get(chrom, 0)
            if pos <= blocked_end:
                if end > blocked_end:
                    missing_positions[chrom] = missing_positions.get(chrom, 0) + callable_bases_in_range(
                        intervals, blocked_end, end)
                    excluded_until[chrom] = end
                skipped += 1
                continue
            if not is_callable(intervals, pos):
                continue
            bases = [rec.REF] + list(rec.ALT)
            if any(len(base) != 1 or base not in "ACGT" for base in bases):
                # Unknown/indel spans are not counted as callable sequence.
                missing_positions[chrom] = missing_positions.get(chrom, 0) + callable_bases_in_range(
                    intervals, pos - 1, end)
                excluded_until[chrom] = end
                skipped += 1
                continue
            choices = []
            for gt in rec.genotypes:
                if len(gt) != 3:
                    raise ValueError(f"Expected diploid genotype at {chrom}:{pos}")
                a, b, phased = gt
                if a is None or b is None or min(a, b) < 0:
                    choices = []
                    break
                if max(a, b) >= len(bases):
                    raise ValueError(f"Invalid allele index at {chrom}:{pos}")
                pair = bases[a] + bases[b]
                choices.append([pair] if phased or a == b else [pair, pair[::-1]])
            if not choices:
                missing_positions[chrom] = missing_positions.get(chrom, 0) + 1
                skipped += 1
                continue
            configurations = sorted({"".join(pairs) for pairs in itertools.product(*choices)})
            if len(set(configurations[0])) < 2:
                continue
            since = callable_bases_in_range(intervals, previous_emitted.get(chrom, 0), pos)
            since -= missing_positions.pop(chrom, 0)
            if since < 1:
                raise ValueError(f"No callable bases before segregating site {chrom}:{pos}")
            previous_emitted[chrom] = pos
            key = chrom if split else None
            if key not in handles:
                path = Path(output) / f"contig_{len(paths):04d}.multihetsep.txt"
                paths.append(str(path))
                handles[key] = stack.enter_context(open(path, "w"))
            handles[key].write(f"{chrom}\t{pos}\t{since}\t{','.join(configurations)}\n")
    if not previous_emitted:
        raise ValueError("No callable segregating SNPs in selected samples")
    print(f"[multihetsep] samples={','.join(vcf.samples)}; contigs={len(previous_emitted)}; "
          f"skipped_uncalled_or_non_SNP={skipped}", file=sys.stderr)
    return paths, list(vcf.samples)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf")
    parser.add_argument("out_multihetsep")
    parser.add_argument("--max-samples", type=int, default=2)
    parser.add_argument("--samples", nargs="+")
    parser.add_argument("--mask")
    parser.add_argument("--assume-callable", action="store_true")
    parser.add_argument("--split", action="store_true")
    args = parser.parse_args()
    if args.max_samples < 1:
        parser.error("--max-samples must be positive")
    selected = args.samples
    if selected is None:
        selected = VCF(args.vcf).samples[:args.max_samples]
    convert(args.vcf, args.out_multihetsep, selected, args.mask, args.assume_callable, args.split)


if __name__ == "__main__":
    main()
