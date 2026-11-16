/*
 * RAiSD
 * Status: stub, explicitly skipped - no conda/PyPI package, and the only
 * source build path depends on the GNU Scientific Library. Per explicit
 * instruction, not worth sinking time into compiling this one; revisit
 * only if a prebuilt binary or container turns up.
 */
process RAISD {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.raisd.stub"), emit: result

    script:
    """
    echo "RAiSD is intentionally skipped: ${meta.reason}" >&2
    exit 1
    """
}
