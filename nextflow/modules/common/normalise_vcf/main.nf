/*
 * Shared per-challenge VCF normalisation (bin/normalise_vcf.sh).
 *
 * The GHIST VCF drops are not clean biallelic call sets: stacked records at
 * one coordinate and triallelic sites occur, and one header line is
 * unterminated. Every VCF-based tool consumes this process's output instead
 * of the raw file, so the per-tool converters only see biallelic SNPs. The
 * excluded positions are emitted as a BED so callable-length-aware tools
 * (PSMC, MSMC2, fastsimcoal2) can subtract them.
 *
 * Only the small BED/summary files are published; the normalised VCF stays
 * in the work directory.
 */
process NORMALISE_VCF {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/normalise_vcf" }, mode: 'copy',
        pattern: '*.{bed,tsv}'
    conda 'bioconda::bcftools bioconda::htslib'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta),
          path("${meta.challenge}.biallelic.vcf.gz"),
          path("${meta.challenge}.biallelic.vcf.gz.tbi"),
          path("${meta.challenge}.excluded_sites.bed"), emit: vcf
    tuple val(meta), path("${meta.challenge}.normalise.tsv"), emit: summary

    script:
    """
    normalise_vcf.sh ${input_file} ${meta.challenge} ${task.cpus}
    """
}
