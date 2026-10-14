/*
 * SHAPEIT5
 * Status: runnable (bioconda: shapeit5)
 *
 * Phasing prerequisite for haplotype-based IBD/relatedness callers (not
 * itself a relatedness estimator, but cataloged alongside them in
 * docs/existingtools.md). No genetic map is provided in this data drop,
 * so this falls back to a constant recombination rate rather than a real
 * map - fine for a survey pass, not for a final submission.
 */
process SHAPEIT5 {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/shapeit5" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.phased.bcf"), emit: result

    script:
    """
    bgzip -c ${input_file} > ${meta.challenge}.vcf.gz
    tabix -p vcf ${meta.challenge}.vcf.gz
    REGION=\$(bcftools view -h ${meta.challenge}.vcf.gz | awk -F'[<,=>]' '/^##contig/{print \$3; exit}')
    phase_common ${meta.args} --input ${meta.challenge}.vcf.gz --region "\$REGION" --output ${meta.challenge}.phased.bcf
    """
}
