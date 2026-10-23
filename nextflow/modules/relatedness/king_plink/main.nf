/*
 * KING + PLINK2
 * Status: runnable (bioconda: plink2)
 *
 * Converts VCF -> PLINK2 pgen/pvar/psam, then runs KING's kinship
 * estimator via plink2's bundled KING-robust algorithm. The challenge
 * supplies the number of diploid autosome pairs (22 human, 38 dog).
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
    plink2 --vcf ${input_file} --make-pgen --out ${meta.challenge} --allow-extra-chr --chr-set ${meta.autosomes} --threads ${task.cpus} --memory ${task.memory.toMega()}
    plink2 --pfile ${meta.challenge} --make-king-table ${meta.args} --out ${meta.challenge} --allow-extra-chr --chr-set ${meta.autosomes} --threads ${task.cpus} --memory ${task.memory.toMega()}
    if [ ! -s ${meta.challenge}.kin0 ]; then
        echo "ERROR: PLINK2 did not produce a nonempty KING table" >&2
        exit 1
    fi
    """
}
