/*
 * NgsRelate
 * Status: runnable (bioconda: ngsrelate, bcftools)
 *
 * NgsRelate natively supports reading genotype calls straight from a VCF
 * (`-T GT -c 1`), bypassing its usual genotype-likelihood route - appropriate
 * here since these are simulated VCFs with hard GT calls, not real
 * low-coverage data with PL/GL fields.
 */
process NGSRELATE {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/ngsrelate" }, mode: 'copy'
    conda { meta.conda }
    // All-pairs EM: ~50k pairs for the 317-sample cohort takes several
    // hours even multi-threaded, well beyond the 2 h default.
    cpus 4
    time '8.h'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.ngsrelate.res"), emit: result

    script:
    """
    ngsRelate -h ${input_file} -T GT -c 1 -p ${task.cpus} -O ${meta.challenge}.ngsrelate.res ${meta.args}
    if [ ! -s ${meta.challenge}.ngsrelate.res ]; then
        echo "ERROR: ngsRelate did not produce a nonempty result" >&2
        exit 1
    fi
    """
}
