/*
 * ADMIXTURE
 * Status: runnable (bioconda: admixture, plink2)
 *
 * ADMIXTURE requires PLINK1 .bed/.bim/.fam with purely numeric chromosome
 * codes; `--allow-extra-chr --chr-set` handles the non-human contigs used
 * by some challenges (e.g. Drosophila 2L/2R/3L/3R, dog autosomes).
 * K is a config-level guess (args: '-K 3 --cv') - tune per challenge/
 * expected population count once real results come back.
 */
process ADMIXTURE {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/admixture" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.*.Q"), emit: result

    script:
    """
    plink2 --vcf ${input_file} --make-bed --out ${meta.challenge} --allow-extra-chr --chr-set -1
    admixture ${meta.args} ${meta.challenge}.bed ${meta.k}
    """
}
