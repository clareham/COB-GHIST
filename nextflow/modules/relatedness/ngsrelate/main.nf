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

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.ngsrelate.res"), emit: result

    script:
    """
    bgzip -c ${input_file} > ${meta.challenge}.vcf.gz
    tabix -p vcf ${meta.challenge}.vcf.gz
    ngsRelate -h ${meta.challenge}.vcf.gz -T GT -c 1 -p ${task.cpus} -O ${meta.challenge}.ngsrelate.res ${meta.args}
    if [ ! -s ${meta.challenge}.ngsrelate.res ]; then
        echo "ERROR: ngsRelate did not produce a nonempty result" >&2
        exit 1
    fi
    """
}
