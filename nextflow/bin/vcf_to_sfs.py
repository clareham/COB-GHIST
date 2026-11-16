#!/usr/bin/env python3
"""Compute a single-population folded site frequency spectrum (SFS) from a
VCF, in fastsimcoal2's *_MAFpop0.obs text format.

SIMPLIFICATION: the GHIST data drop doesn't include a population map for
the two-population challenges (secondary_contact, split_migration), so
this script pools every sample into one population and emits a 1D folded
SFS rather than fastsimcoal2's 2D joint SFS format. That's enough to run
fastsimcoal2 end to end as a single-population model for an early survey,
but a faithful 2-population analysis needs a popmap (sample -> pop1/pop2)
added once that information is available, feeding a joint-SFS variant of
this script.
"""
import argparse

from cyvcf2 import VCF


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf")
    ap.add_argument("out_obs")
    args = ap.parse_args()

    vcf = VCF(args.vcf)
    n_hap = 2 * len(vcf.samples)
    n_fold = n_hap // 2  # folded spectrum has floor(n_hap/2) + 1 bins (0..n_fold)
    counts = [0] * (n_fold + 1)

    for rec in vcf:
        ac = 0
        an = 0
        for gt in rec.genotypes:
            for a in gt[:2]:
                if a is None or a < 0:
                    continue
                an += 1
                if a != 0:
                    ac += 1
        if an == 0:
            continue
        maf = min(ac, an - ac)
        bin_idx = min(maf, n_fold)
        counts[bin_idx] += 1

    with open(args.out_obs, "w") as out:
        out.write("1 observations\n")
        out.write("\t".join(f"d0_{i}" for i in range(n_fold + 1)) + "\n")
        out.write("\t".join(str(c) for c in counts) + "\n")


if __name__ == "__main__":
    main()
