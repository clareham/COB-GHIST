#!/usr/bin/env python3
"""Build an MSMC2 'multihetsep' file directly from a VCF.

msmc-tools' own generate_multihetsep.py isn't packaged on bioconda/PyPI
(it's a small script you'd otherwise have to git-clone by hand), so this
reproduces the multihetsep format directly:

  --mask not given: every site is treated as callable and
    'nSitesSincePrevCall' is just the physical distance since the
    previous emitted site. Appropriate for the regular demography VCF
    challenges, where no reference/coverage mask exists at all.

  --mask <bed>: 'nSitesSincePrevCall' instead counts actual callable
    bases (from a real callability mask, e.g. bin/bam_to_depth_mask.py's
    output for the low-coverage BAM-based challenges), and sites falling
    outside any masked interval are skipped entirely - this is the
    correct, mask-aware behaviour MSMC2 expects.

Output columns (tab-separated): CHROM POS N_SITES_SINCE_PREV ALLELES
"""
import argparse
import bisect
from collections import defaultdict

from cyvcf2 import VCF


def load_mask(path):
    """BED (0-based, half-open) -> {chrom: (starts, ends)}, sorted & assumed non-overlapping."""
    by_chrom = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split()[:3]
            by_chrom[chrom].append((int(start), int(end)))
    out = {}
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        out[chrom] = ([s for s, _ in ivs], [e for _, e in ivs])
    return out


def callable_bases_in_range(mask_chrom, start_excl, end_incl):
    """Count callable (masked-in) bases in the 1-based range (start_excl, end_incl]."""
    if mask_chrom is None:
        return end_incl - start_excl
    starts, ends = mask_chrom
    total = 0
    # intervals are 0-based [start, end); a 1-based position p is covered if start < p <= end
    i = bisect.bisect_right(ends, start_excl)
    while i < len(starts) and starts[i] < end_incl:
        total += min(ends[i], end_incl) - max(starts[i], start_excl)
        i += 1
    return max(total, 0)


def is_callable(mask_chrom, pos_1based):
    if mask_chrom is None:
        return True
    starts, ends = mask_chrom
    i = bisect.bisect_right(starts, pos_1based - 1) - 1
    return 0 <= i < len(starts) and starts[i] < pos_1based <= ends[i]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf")
    ap.add_argument("out_multihetsep")
    ap.add_argument("--max-samples", type=int, default=2,
                     help="Number of diploid samples (individuals) to include (default: 2, i.e. 4 haplotypes)")
    ap.add_argument("--mask", default=None,
                     help="Callability mask BED - restricts callable sites and makes nSitesSincePrevCall exact")
    args = ap.parse_args()

    mask = load_mask(args.mask) if args.mask else None

    vcf = VCF(args.vcf)
    samples = vcf.samples[:args.max_samples]
    vcf = VCF(args.vcf, samples=samples)

    with open(args.out_multihetsep, "w") as out:
        prev_pos = {}
        for rec in vcf:
            mask_chrom = mask.get(rec.CHROM) if mask is not None else None
            if mask is not None and not is_callable(mask_chrom, rec.POS):
                continue
            alleles = []
            any_variant = False
            ref, alts = rec.REF, rec.ALT
            bases = [ref] + list(alts)
            ok = True
            for gt in rec.genotypes:
                a, b = gt[0], gt[1]
                if a is None or b is None or a < 0 or b < 0 or a >= len(bases) or b >= len(bases):
                    ok = False
                    break
                if a != 0 or b != 0:
                    any_variant = True
                alleles.append(bases[a])
                alleles.append(bases[b])
            if not ok or not any_variant:
                continue
            if any(len(x) != 1 for x in alleles):
                continue  # skip multi-base indels, multihetsep wants single bases
            prev = prev_pos.get(rec.CHROM, 0)
            n_since = callable_bases_in_range(mask_chrom, prev, rec.POS) if mask is not None else rec.POS - prev
            prev_pos[rec.CHROM] = rec.POS
            out.write(f"{rec.CHROM}\t{rec.POS}\t{n_since}\t{''.join(alleles)}\n")


if __name__ == "__main__":
    main()
