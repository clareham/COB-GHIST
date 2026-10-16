#!/usr/bin/env bash
# Joint-call two diploid samples and intersect their per-base callability.
# Usage: <archive.tar.gz> <ref.fasta> <out_prefix> [min_depth] [max_depth]
set -euo pipefail

ARCHIVE=$1
REF=$2
OUT_PREFIX=$3
MIN_DEPTH=${4:-2}
MAX_DEPTH=${5:-20}
BIN_DIR=$(CDPATH= cd -- "$(dirname -- "$(readlink -f "$0")")" && pwd)

samtools faidx "$REF"
python3 "$BIN_DIR/bam_archive_samples.py" "$ARCHIVE" extracted "${OUT_PREFIX}.populations.tsv"
mapfile -t BAMS < <(tail -n +2 "${OUT_PREFIX}.populations.tsv" | cut -f3)
if [ "${#BAMS[@]}" -ne 2 ]; then
    echo "ERROR: expected two selected BAM samples" >&2
    exit 1
fi

MASKS=()
for index in 0 1; do
    bam=${BAMS[$index]}
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam%.bam}.bai" ]; then
        samtools index "$bam"
    fi
    mask="${OUT_PREFIX}.sample${index}.mask.bed"
    python3 "$BIN_DIR/bam_to_depth_mask.py" "$bam" "$mask" --min-depth "$MIN_DEPTH" --max-depth "$MAX_DEPTH"
    MASKS+=("$mask")
done

bcftools mpileup -q 20 -Q 20 -C50 -Ou -f "$REF" "${BAMS[@]}" \
    | bcftools call -m -v -V indels -Oz -o "${OUT_PREFIX}.merged.vcf.gz"
tabix -p vcf "${OUT_PREFIX}.merged.vcf.gz"
bedtools intersect -a "${MASKS[0]}" -b "${MASKS[1]}" \
    | bedtools sort -i - | bedtools merge -i - > "${OUT_PREFIX}.mask.bed"
if [ ! -s "${OUT_PREFIX}.mask.bed" ]; then
    echo "ERROR: selected individuals have no jointly callable bases" >&2
    exit 1
fi
