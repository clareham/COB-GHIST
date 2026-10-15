/*
 * Flex-Sweep - disabled integration (bioconda package is available).
 *
 * Requires feature extraction and a trained classifier connected through
 * Flex-Sweep's subcommands, not a top-level --vcf/--output invocation.
 * See https://flexsweep.readthedocs.io/en/latest.
 */
process FLEXSWEEP {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.flexsweep.out"), emit: result

    script:
    """
    echo "ERROR: Flex-Sweep integration is disabled: implement feature extraction and trained-model prediction before enabling it." >&2
    exit 1
    """
}
