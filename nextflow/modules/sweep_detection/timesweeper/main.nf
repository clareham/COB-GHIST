/*
 * Timesweeper
 * Status: excluded, not a packaging problem - Timesweeper is pip-
 * installable, but it's designed for serial/time-series sampling (e.g.
 * ancient + modern individuals), and every GHIST sweep-detection testing
 * VCF here is a single present-day sampling (see onesweep.config etc. -
 * sample IDs are generic tsk_N with no time metadata). There's no time
 * axis for it to use on this data as provided.
 */
process TIMESWEEPER {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.timesweeper.stub"), emit: result

    script:
    """
    echo "Timesweeper is excluded for this data: ${meta.reason}" >&2
    exit 1
    """
}
