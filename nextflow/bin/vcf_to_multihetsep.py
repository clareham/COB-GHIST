#!/usr/bin/env python3
"""Build an MSMC2 'multihetsep' file directly from a VCF.

msmc-tools' own generate_multihetsep.py isn't packaged on bioconda/PyPI
(it's a small script you'd otherwise have to git-clone by hand), and - as
with the PSMC path - there's no reference/mask BAM/coverage info in this
data drop to compute a proper callability mask. This script reproduces
the multihetsep format for the common case here: unphased simulated
genotypes with no missingness, so every site is treated as callable and
'nSitesSincePrevCall' is just the physical distance since the previous
emitted site.

CAVEAT: without a real strict/callability mask, MSMC2's effective
population size scaling will be off versus a production run - fine for an
early tool survey, not for the final submission.

Output columns (tab-separated): CHROM POS N_SITES_SINCE_PREV ALLELES
"""
import argparse

from cyvcf2 import VCF


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf")
    ap.add_argument("out_multihetsep")
    ap.add_argument("--max-samples", type=int, default=2,
                     help="Number of diploid samples (individuals) to include (default: 2, i.e. 4 haplotypes)")
    args = ap.parse_args()

    vcf = VCF(args.vcf)
    samples = vcf.samples[:args.max_samples]
    vcf = VCF(args.vcf, samples=samples)

    with open(args.out_multihetsep, "w") as out:
        prev_pos = {}
        for rec in vcf:
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
            n_since = rec.POS - prev
            prev_pos[rec.CHROM] = rec.POS
            out.write(f"{rec.CHROM}\t{rec.POS}\t{n_since}\t{''.join(alleles)}\n")


if __name__ == "__main__":
    main()
