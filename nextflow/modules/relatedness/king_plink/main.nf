/*
 * KING + PLINK2
 * Status: runnable (bioconda: king, plink2)
 *
 * Converts VCF -> PLINK2 pgen/pvar/psam, then runs KING's kinship
 * estimator (`--kinship`) via plink2's bundled KING-robust algorithm, so
 * only one binary round-trip is needed (no intermediate PLINK1 .bed).
 */
process KING_PLINK {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/king_plink" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.kin0"), emit: result

    script:
    """
    plink2 --vcf ${input_file} --make-pgen --out ${meta.challenge} --allow-extra-chr --chr-set -1
    plink2 --pfile ${meta.challenge} --make-king-table ${meta.args} --out ${meta.challenge} --allow-extra-chr --chr-set -1
    """
}
