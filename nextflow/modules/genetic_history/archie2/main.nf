/*
 * ArchIE2
 * Status: stub - introgressed-DNA detection tool, tangential fit for the
 * demography challenges (noted as "may be tangentially useful" in
 * docs/existingtools.md). No conda/PyPI package found. To activate:
 * compile from source manually, flip `enabled: true` deliberately (it's
 * disabled even once installed, pending a clearer fit to secondary_contact),
 * and fill in the real invocation.
 */
process ARCHIE2 {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.archie2.stub"), emit: result

    script:
    """
    echo "ArchIE2 is not wired up yet: ${meta.reason}" >&2
    exit 1
    """
}
