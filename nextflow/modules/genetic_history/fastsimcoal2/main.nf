/*
 * fastsimcoal2 - disabled integration (bioconda package is available).
 *
 * Requires challenge-specific .tpl/.est models and a verified folded-SFS
 * invocation with sample sizes matching the observed spectrum. Two-
 * population challenges also need a population map and joint SFS.
 * The former constant-size, one-sample placeholder was not a valid
 * implementation of these challenges.
 */
process FASTSIMCOAL2 {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.fsc2_out"), emit: result

    script:
    """
    echo "ERROR: fastsimcoal2 integration is disabled: implement challenge-specific models and matching SFS input before enabling it." >&2
    exit 1
    """
}
