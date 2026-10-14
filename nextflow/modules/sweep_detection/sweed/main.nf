/*
 * SweeD
 * Status: stub - no conda/PyPI package found, needs a from-source C build.
 * To activate: compile manually, flip `enabled: true`, and fill in the
 * real invocation (SweeD takes a standard VCF or ms-style input directly,
 * simpler than SweepFinder2's custom .freq format).
 */
process SWEED {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.sweed.stub"), emit: result

    script:
    """
    echo "SweeD is not wired up yet: ${meta.reason}" >&2
    exit 1
    """
}
