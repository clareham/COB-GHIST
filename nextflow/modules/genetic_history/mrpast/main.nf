/*
 * mrpast - disabled integration (bioconda package is available).
 *
 * Requires a demographic model, sample-to-population mapping, dated ARGs,
 * and the documented `mrpast process` / `mrpast solve` workflow.
 * Raw tsinfer output alone is not a substitute for dated ARG inference.
 * See https://mrpast.readthedocs.io/en/latest/workflows/tsinfer.html.
 */
process MRPAST {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.mrpast_out"), emit: result

    script:
    """
    echo "ERROR: mrpast integration is disabled: implement model/population mapping, dated ARG inference, and process/solve before enabling it." >&2
    exit 1
    """
}
