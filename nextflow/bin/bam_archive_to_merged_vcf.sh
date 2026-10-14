#!/usr/bin/env bash
# Build a merged multi-sample VCF + intersected callability mask from a
# GHIST low-coverage BAM archive, for MSMC2's multihetsep input
# (bin/vcf_to_multihetsep.py --mask). Picks up to 2 individuals: one per
# population if the archive's popfile.txt lists >1 population (giving
# MSMC2 the cross-population haplotype pairing it needs for a split-time
# estimate), otherwise the first 2 individuals from the single population.
#
# Usage: bam_archive_to_merged_vcf.sh <archive.tar.gz> <ref.fasta> <out_prefix>
# Outputs: <out_prefix>.merged.vcf.gz (+ .tbi), <out_prefix>.mask.bed
set -euo pipefail

ARCHIVE=$1
REF=$2
OUT_PREFIX=$3

samtools faidx "$REF"

mkdir -p extracted
tar xzf "$ARCHIVE" -C extracted
ARCHIVE_DIR=$(find extracted -mindepth 1 -maxdepth 1 -type d | head -1)
POPFILE="${ARCHIVE_DIR}/popfile.txt"

# Pick up to 2 samples: prefer one per distinct population.
mapfile -t SAMPLES < <(python3 - "$POPFILE" <<'PY'
import sys
# popfile.txt has a header row ('sample\tpop') - skip it, since otherwise
# the literal header text gets treated as a (nonexistent) sample/pop pair.
rows = [l.split() for l in open(sys.argv[1]).readlines()[1:] if l.strip()]
pops = {}
for sample, pop in rows:
    pops.setdefault(pop, sample)
picks = list(pops.values())[:2]
if len(picks) < 2:
    picks = [r[0] for r in rows[:2]]
print("\n".join(picks))
PY
)

find_bam() {
    local sample=$1
    for f in "${ARCHIVE_DIR}"/*.bam; do
        sm=$(samtools view -H "$f" | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if ($i ~ /^SM:/) print substr($i,4)}' | head -1)
        if [ "$sm" = "$sample" ]; then echo "$f"; return 0; fi
    done
    return 1
}

VCFS=()
MASKS=()
for sample in "${SAMPLES[@]}"; do
    bam=$(find_bam "$sample")
    if [ -z "$bam" ]; then
        echo "ERROR: no BAM with read-group SM matching sample '$sample'" >&2
        exit 1
    fi
    bcftools mpileup -C50 -f "$REF" "$bam" | bcftools call -c -V indels -Oz -o "${sample}.vcf.gz"
    tabix -p vcf "${sample}.vcf.gz"
    bam_to_depth_mask.py "$bam" "${sample}.mask.bed"
    VCFS+=("${sample}.vcf.gz")
    MASKS+=("${sample}.mask.bed")
done

bcftools merge "${VCFS[@]}" -Oz -o "${OUT_PREFIX}.merged.vcf.gz"
tabix -p vcf "${OUT_PREFIX}.merged.vcf.gz"

if [ "${#MASKS[@]}" -gt 1 ]; then
    bedtools intersect -a "${MASKS[0]}" -b "${MASKS[@]:1}" > "${OUT_PREFIX}.mask.bed"
else
    cp "${MASKS[0]}" "${OUT_PREFIX}.mask.bed"
fi
