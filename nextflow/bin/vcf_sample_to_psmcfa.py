#!/usr/bin/env python3
"""Build PSMC input from a diploid variant-only VCF.

Requires a callable BED or explicit --assume-callable for simulated
error-free genomes. Header lengths retain invariant and trailing sequence.
Only fully callable windows are retained; missing genotypes mask their
whole window conservatively. T=homozygous, K=heterozygous, N=uncallable.
"""

import argparse
import sys

from cyvcf2 import VCF

from callable_regions import callable_bases_in_range, contig_lengths, load_mask


def convert(vcf_path, output, sample=None, window=100, mask_path=None, assume_callable=False):
    if window <= 0:
        raise ValueError("Window size must be positive")
    if mask_path is None and not assume_callable:
        raise ValueError("Provide --mask or explicitly choose --assume-callable")
    vcf = VCF(vcf_path)
    if not vcf.samples:
        raise ValueError("VCF has no samples")
    sample = sample or vcf.samples[0]
    if sample not in vcf.samples:
        raise ValueError(f"Sample not in VCF: {sample}")
    vcf.set_samples([sample])
    lengths = contig_lengths(vcf)
    mask = load_mask(mask_path) if mask_path else None
    sequences = {}
    missing = {}
    for chrom, length in lengths.items():
        sequence = bytearray(b"N" * ((length + window - 1) // window))
        intervals = mask.get(chrom, ([], [])) if mask is not None else None
        for index in range(len(sequence)):
            start, end = index * window, min((index + 1) * window, length)
            if callable_bases_in_range(intervals, start, end) == end - start:
                sequence[index] = ord("T")
        sequences[chrom] = sequence
        missing[chrom] = set()
    previous = {}
    for rec in vcf:
        if rec.CHROM not in lengths or not 1 <= rec.POS <= lengths[rec.CHROM]:
            raise ValueError(f"Record outside declared contig: {rec.CHROM}:{rec.POS}")
        if rec.POS <= previous.get(rec.CHROM, 0):
            raise ValueError("VCF positions must be strictly increasing within each contig")
        previous[rec.CHROM] = rec.POS
        genotype = rec.genotypes[0]
        index = (rec.POS - 1) // window
        if len(genotype) != 3:
            raise ValueError(f"Expected diploid genotype at {rec.CHROM}:{rec.POS}")
        alleles = genotype[:2]
        if any(allele is None or allele < 0 for allele in alleles):
            missing[rec.CHROM].add(index)
            continue
        if any(allele >= 1 + len(rec.ALT) for allele in alleles):
            raise ValueError(f"Invalid genotype at {rec.CHROM}:{rec.POS}")
        if (len(rec.REF) != 1 or rec.REF not in "ACGT" or
                any(len(alt) != 1 or alt not in "ACGT" for alt in rec.ALT)):
            end = min(rec.POS + len(rec.REF) - 1, lengths[rec.CHROM])
            for affected in range(index, (end - 1) // window + 1):
                missing[rec.CHROM].add(affected)
            continue
        if alleles[0] != alleles[1] and sequences[rec.CHROM][index] != ord("N"):
            sequences[rec.CHROM][index] = ord("K")
    heterozygous = 0
    for chrom, sequence in sequences.items():
        for index in missing[chrom]:
            sequence[index] = ord("N")
        heterozygous += sequence.count(b"K")
    if heterozygous == 0:
        raise ValueError(f"Sample {sample} has no callable heterozygous PSMC windows")
    with open(output, "w") as handle:
        for chrom, sequence in sequences.items():
            if all(base == ord("N") for base in sequence):
                continue
            handle.write(f">{chrom}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index:index + 60].decode("ascii") + "\n")
    print(f"[psmc input] sample={sample}; heterozygous_windows={heterozygous}; "
          f"callability={'BED' if mask_path else 'assumed across declared contigs'}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf")
    parser.add_argument("out_psmcfa")
    parser.add_argument("--sample")
    parser.add_argument("--window", type=int, default=100)
    parser.add_argument("--mask")
    parser.add_argument("--assume-callable", action="store_true")
    args = parser.parse_args()
    convert(args.vcf, args.out_psmcfa, args.sample, args.window, args.mask, args.assume_callable)


if __name__ == "__main__":
    main()
