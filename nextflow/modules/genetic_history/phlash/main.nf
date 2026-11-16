/*
 * PHLASH (Population History Learning by Averaging Sampled Histories)
 * Status: stub - Python package installed from the GitHub source tree,
 * no PyPI/conda package found. To activate: `pip install` from the
 * project's GitHub repo inside a conda env, point `conda`/`container` at
 * that environment, flip `enabled: true`, and fill in the real PHLASH
 * Python-API call below (it works off tree sequences / ARGs, similar to
 * mrpast - see modules/genetic_history/mrpast/main.nf for the tsinfer
 * conversion pattern).
 */
process PHLASH {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.phlash.stub"), emit: result

    script:
    """
    echo "PHLASH is not wired up yet: ${meta.reason}" >&2
    exit 1
    """
}
