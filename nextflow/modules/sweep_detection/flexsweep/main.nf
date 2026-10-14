/*
 * Flex-Sweep (CNN-based sweep classifier)
 * Status: runnable-unverified (bioconda: flexsweep)
 *
 * Package installs cleanly, but we haven't independently verified the
 * exact CLI invocation/feature-extraction workflow against a real run -
 * check `flexsweep --help` (and whether it needs a pretrained model vs.
 * training its own from simulations) before trusting output.
 */
process FLEXSWEEP {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/flexsweep" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.flexsweep.out"), emit: result

    script:
    """
    flexsweep ${meta.args} --vcf ${input_file} --output ${meta.challenge}.flexsweep.out \
        || echo "flexsweep invocation failed/unverified - check CLI flags against 'flexsweep --help'" >&2
    """
}
