/*
 * IcMLkin
 * Status: stub - low-coverage maximum-likelihood kinship estimator, no
 * conda/PyPI package found. To activate: compile from source manually,
 * flip `enabled: true`, and fill in the real invocation (needs raw reads,
 * not provided for these VCF-based relatedness challenges - best suited
 * to a future low-coverage relatedness challenge if one is added).
 */
process ICMLKIN {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.icmlkin.stub"), emit: result

    script:
    """
    echo "IcMLkin is not wired up yet: ${meta.reason}" >&2
    exit 1
    """
}
