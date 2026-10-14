#!/usr/bin/env python3
"""Convert a VCF to SweepFinder2's .freq input format.

Columns: position  x  n  folded
  x      = derived allele count at the site
  n      = number of called haplotypes (chromosomes) at the site
  folded = 0 (unfolded) or 1 (folded)

ASSUMPTION: these are tskit/SLiM simulated VCFs, where REF is written as
the true ancestral allele by convention - so we treat ALT count as the
derived allele count and emit an *unfolded* spectrum (folded=0). If that
assumption doesn't hold for a given scenario, rerun with --folded to fall
back to minor-allele counts instead.
"""
import argparse

from cyvcf2 import VCF


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf")
    ap.add_argument("out_freq")
    ap.add_argument("--folded", action="store_true", help="Emit folded (minor-allele) spectrum instead of unfolded")
    args = ap.parse_args()

    with open(args.out_freq, "w") as out:
        out.write("position\tx\tn\tfolded\n")
        for rec in VCF(args.vcf):
            ac = 0
            an = 0
            for gt in rec.genotypes:
                for a in gt[:2]:
                    if a is None or a < 0:
                        continue
                    an += 1
                    if a != 0:
                        ac += 1
            if an == 0 or ac == 0 or ac == an:
                continue  # SweepFinder2 skips invariant sites
            x = min(ac, an - ac) if args.folded else ac
            out.write(f"{rec.POS}\t{x}\t{an}\t{1 if args.folded else 0}\n")


if __name__ == "__main__":
    main()
