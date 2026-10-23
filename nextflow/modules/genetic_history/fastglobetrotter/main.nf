/*
 * fastGLOBETROTTER
 * Status: stub - linkage-disequilibrium-based admixture dating method, no
 * conda/PyPI package found. To activate: compile from source manually,
 * flip `enabled: true` in secondary_contact/split_migration configs, and
 * fill in the real invocation (needs a painted-chromosome / ChromoPainter
 * output as input, which itself needs a separate ChromoPainter run first).
 */
process FASTGLOBETROTTER {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.fastglobetrotter.stub"), emit: result

    script:
    """
    echo "fastGLOBETROTTER is not wired up yet: ${meta.reason}" >&2
    exit 1
    """
}
