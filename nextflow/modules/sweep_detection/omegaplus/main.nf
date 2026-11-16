/*
 * OmegaPlus
 * Status: stub - no conda/PyPI package found, needs a from-source C build.
 * To activate: compile manually, flip `enabled: true`, and fill in the
 * real invocation (OmegaPlus takes a VCF/ms-style alignment directly and
 * computes a composite omega-statistic across sliding windows).
 */
process OMEGAPLUS {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.omegaplus.stub"), emit: result

    script:
    """
    echo "OmegaPlus is not wired up yet: ${meta.reason}" >&2
    exit 1
    """
}
