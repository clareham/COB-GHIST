/*
 * MSMC2 (Multiple Sequentially Markovian Coalescent)
 * Status: runnable (bioconda: msmc2, cyvcf2)
 *
 * Builds a multihetsep file directly from genotypes
 * (bin/vcf_to_multihetsep.py) since msmc-tools' own generator isn't
 * packaged and there's no coverage/callability mask in this data drop -
 * see that script's docstring for the exact caveat (every site treated
 * as callable).
 *
 * Not offered for the *_low_coverage challenges (see their challenge
 * configs) - those hand out a single haploid sequence with no het
 * signal, which MSMC2 cannot use at all.
 */
process MSMC2 {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/msmc2" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.msmc2.final.txt"), emit: result

    script:
    """
    vcf_to_multihetsep.py ${input_file} ${meta.challenge}.multihetsep.txt
    msmc2 ${meta.args} -o ${meta.challenge}.msmc2 ${meta.challenge}.multihetsep.txt
    """
}
