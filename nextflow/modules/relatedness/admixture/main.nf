/*
 * ADMIXTURE
 * Status: runnable (bioconda: admixture, plink2)
 *
 * ADMIXTURE requires PLINK1 .bed/.bim/.fam. The challenge supplies the
 * diploid autosome count for PLINK (22 human, 38 dog).
 * K is a config-level guess (default: 3) - tune per challenge/
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
    plink2 --vcf ${input_file} --make-bed --out ${meta.challenge} --allow-extra-chr --chr-set ${meta.autosomes} --threads ${task.cpus} --memory ${task.memory.toMega()}
    admixture -j${task.cpus} ${meta.args} ${meta.challenge}.bed ${meta.k}
    if [ ! -s ${meta.challenge}.${meta.k}.Q ]; then
        echo "ERROR: ADMIXTURE did not produce a nonempty ancestry matrix" >&2
        exit 1
    fi
    """
}
