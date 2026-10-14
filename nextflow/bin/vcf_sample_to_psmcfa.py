#!/usr/bin/env python3
"""Build a PSMC-fa file for one sample directly from a VCF.

The GHIST simulated challenges ship no reference FASTA for most
demographic-history scenarios (the contigs are purely simulated), so the
usual `bcftools consensus -f ref.fa | vcfutils.pl vcf2fq | fq2psmcfa`
route doesn't apply. Classic PSMC only needs, per non-overlapping window,
whether the chosen sample is heterozygous anywhere in that window - so we
derive the 'T'/'K'/'N' call string straight from genotypes.

Encoding (matches fq2psmcfa's own output alphabet):
  T = window has >=1 called site, all homozygous
  K = window has >=1 called site, at least one heterozygous
  N = window has no called sites at all
"""
import argparse
import sys

from cyvcf2 import VCF


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf")
    ap.add_argument("out_psmcfa")
    ap.add_argument("--sample", default=None, help="Sample ID (default: first sample in the VCF)")
    ap.add_argument("--window", type=int, default=100, help="Window size in bp (PSMC default: 100)")
    args = ap.parse_args()

    vcf = VCF(args.vcf, samples=[args.sample] if args.sample else None)
    sample = vcf.samples[0]
    print(f"[vcf_sample_to_psmcfa] using sample '{sample}'", file=sys.stderr)

    with open(args.out_psmcfa, "w") as out:
        contig = None
        bins = {}
        max_bin = -1

        def flush(contig_name):
            if contig_name is None:
                return
            out.write(f">{contig_name}\n")
            seq = "".join(bins.get(i, "N") for i in range(max_bin + 1))
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")

        for rec in vcf:
            if contig != rec.CHROM:
                flush(contig)
                contig, bins, max_bin = rec.CHROM, {}, -1
            gt = rec.genotypes[0]
            alleles = [a for a in gt[:2] if a is not None and a >= 0]
            if len(alleles) < 2:
                continue  # uncalled site, contributes nothing (stays 'N' unless another site fills the bin)
            b = (rec.POS - 1) // args.window
            max_bin = max(max_bin, b)
            het = alleles[0] != alleles[1]
            if bins.get(b) == "K":
                continue
            bins[b] = "K" if het else "T"
        flush(contig)


if __name__ == "__main__":
    main()
