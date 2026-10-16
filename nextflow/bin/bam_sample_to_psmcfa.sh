#!/usr/bin/env bash
# Classic PSMC input-prep pipeline for one low-coverage BAM sample, per
# PSMC's own documented usage (lh3/psmc README) - this is the *real*
# intended PSMC input: low-coverage reads aligned to a reference, called
# into a per-site diploid consensus that encodes heterozygosity via
# IUPAC ambiguity codes, which PSMC then reads directly.
#
# Unlike bin/vcf_sample_to_psmcfa.py (used when no reference genome
# exists at all, e.g. the regular demography VCF challenges), this script
# needs a real reference FASTA + aligned BAM, both of which the GHIST
# *_low_coverage challenges provide (once their BAM archives are
# downloaded - see conf/genetic_history/bottleneck_low_coverage.config).
#
# Usage: bam_sample_to_psmcfa.sh <ref.fasta> <sample.bam> <out_prefix> [min_depth] [max_depth]
set -euo pipefail

REF=$1
BAM=$2
OUT_PREFIX=$3
MIN_DEPTH=${4:-1}
MAX_DEPTH=${5:-10}

samtools faidx "$REF"
bcftools mpileup -q20 -Q20 -C50 -f "$REF" "$BAM" \
    | bcftools call -c \
    | vcfutils.pl vcf2fq -d "$MIN_DEPTH" -D "$MAX_DEPTH" \
    > "${OUT_PREFIX}.consensus.fq"
fq2psmcfa -q20 "${OUT_PREFIX}.consensus.fq" > "${OUT_PREFIX}.psmcfa"
