#!/usr/bin/env python3
"""Convert single-contig biallelic SNPs into SweepFinder2 frequency input.

Use --folded when ancestral alleles are unknown. Unfolded mode requires
valid INFO/AA values; REF is not assumed ancestral.
"""

import argparse
import sys

from cyvcf2 import VCF

from callable_regions import is_callable, iter_positions, load_mask


def convert(vcf_path, output, folded=False, mask_path=None):
    mask = load_mask(mask_path) if mask_path else None
    contig, emitted, skipped = None, 0, 0
    with open(output, "w") as handle:
        handle.write("position\tx\tn\tfolded\n")
        for records in iter_positions(VCF(vcf_path)):
            rec = records[0]
            if contig is None:
                contig = rec.CHROM
            if rec.CHROM != contig:
                raise ValueError("SweepFinder2 requires a separate scan per contig")
            intervals = mask.get(contig, ([], [])) if mask is not None else None
            if not is_callable(intervals, rec.POS):
                skipped += 1
                continue
            if (len(records) > 1 or len(rec.ALT) != 1 or len(rec.REF) != 1 or len(rec.ALT[0]) != 1 or
                    rec.REF not in "ACGT" or rec.ALT[0] not in "ACGT"):
                skipped += 1
                continue
            alternate, count = 0, 0
            for gt in rec.genotypes:
                if len(gt) != 3:
                    raise ValueError("SweepFinder2 survey expects diploid genotypes")
                for allele in gt[:2]:
                    if allele is None or allele < 0:
                        continue
                    if allele > 1:
                        raise ValueError(f"Invalid allele at {rec.CHROM}:{rec.POS}")
                    count += 1
                    alternate += allele
            if count == 0 or alternate in (0, count):
                continue
            if folded:
                derived = min(alternate, count - alternate)
            else:
                ancestor = rec.INFO.get("AA")
                if ancestor not in (rec.REF, rec.ALT[0]):
                    raise ValueError(f"Missing/invalid INFO/AA at {contig}:{rec.POS}; use --folded")
                derived = alternate if ancestor == rec.REF else count - alternate
            handle.write(f"{rec.POS}\t{derived}\t{count}\t{int(folded)}\n")
            emitted += 1
    if emitted == 0:
        raise ValueError("No polymorphic callable SNPs for SweepFinder2")
    print(f"[SweepFinder2 input] contig={contig}; sites={emitted}; "
          f"masked_or_non_biallelic_SNP={skipped}; folded={folded}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf")
    parser.add_argument("out_freq")
    parser.add_argument("--folded", action="store_true")
    parser.add_argument("--mask")
    args = parser.parse_args()
    convert(args.vcf, args.out_freq, args.folded, args.mask)


if __name__ == "__main__":
    main()
