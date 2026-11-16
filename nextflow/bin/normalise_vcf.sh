#!/usr/bin/env bash
# Normalise a GHIST challenge VCF into a clean biallelic-SNP call set.
#
# The simulated drops are not clean call sets: several distinct simulated
# sites can collapse onto one coordinate (stacked records, often with
# different REF bases), triallelic records occur, and the `##INFO=<ID=AA`
# header line is unterminated (htslib otherwise swallows the following
# header lines). This script repairs the header, then removes every position
# that is not exactly one biallelic SNP record. The removed positions are
# written to <prefix>.excluded_sites.bed so callable-length-aware tools can
# subtract them.
#
# Usage: normalise_vcf.sh <in.vcf[.gz]> <out_prefix> [threads]
# Outputs: <prefix>.biallelic.vcf.gz(.tbi), <prefix>.excluded_sites.bed,
#          <prefix>.normalise.tsv
set -euo pipefail

IN=$1
PREFIX=$2
THREADS=${3:-1}

read_vcf() {
    case "$1" in
        *.gz) gzip -cd "$1" ;;
        *) cat "$1" ;;
    esac
}

# 1. Repair unterminated `##KEY=<...>` header lines and convert to BCF once,
#    so the following passes read a compact indexed-free binary stream.
read_vcf "$IN" \
    | awk -v counts="${PREFIX}.repaired.txt" '
        /^##/ && /=</ && !/>$/ {
            quotes = gsub(/"/, "\"", $0)
            $0 = $0 (quotes % 2 ? "\">" : ">")
            repaired++
        }
        { print }
        END { print repaired+0 > counts }' \
    | bcftools view --threads "$THREADS" -Ob -o "${PREFIX}.raw.bcf" -

# 2. Positions to drop: any coordinate carrying more than one record, and any
#    record that is not a biallelic SNP. BED is 0-based half-open.
{
    bcftools query -f '%CHROM\t%POS\n' "${PREFIX}.raw.bcf" \
        | uniq -d | awk -F'\t' '{print $1"\t"$2-1"\t"$2}'
    bcftools view -i 'N_ALT!=1 || TYPE!="snp"' "${PREFIX}.raw.bcf" \
        | bcftools query -f '%CHROM\t%POS0\t%END\n'
} | sort -k1,1 -k2,2n -u > "${PREFIX}.excluded_sites.bed"

# 3. Emit the retained biallelic SNPs. `-T ^file` streams without an index.
if [ -s "${PREFIX}.excluded_sites.bed" ]; then
    bcftools view --threads "$THREADS" -T "^${PREFIX}.excluded_sites.bed" -Oz \
        -o "${PREFIX}.biallelic.vcf.gz" "${PREFIX}.raw.bcf"
else
    bcftools view --threads "$THREADS" -Oz -o "${PREFIX}.biallelic.vcf.gz" "${PREFIX}.raw.bcf"
fi
tabix -p vcf "${PREFIX}.biallelic.vcf.gz"

raw=$(bcftools view -H "${PREFIX}.raw.bcf" | wc -l)
retained=$(bcftools index -n "${PREFIX}.biallelic.vcf.gz")
excluded=$(wc -l < "${PREFIX}.excluded_sites.bed")
repaired=$(cat "${PREFIX}.repaired.txt")
if [ "$retained" -eq 0 ]; then
    echo "ERROR: no biallelic SNPs remain after normalising $IN" >&2
    exit 1
fi
{
    printf 'metric\tvalue\n'
    printf 'raw_records\t%s\n' "$raw"
    printf 'excluded_positions\t%s\n' "$excluded"
    printf 'excluded_records\t%s\n' "$((raw - retained))"
    printf 'retained_biallelic_snps\t%s\n' "$retained"
    printf 'repaired_header_lines\t%s\n' "$repaired"
} > "${PREFIX}.normalise.tsv"
rm -f "${PREFIX}.raw.bcf" "${PREFIX}.repaired.txt"
cat "${PREFIX}.normalise.tsv" >&2
