/*
 * READ + KIN
 * Status: stub - both READ and KIN are standalone GitHub Python scripts,
 * not packaged on conda/PyPI. To activate: clone and vendor them (e.g.
 * under bin/), flip `enabled: true`, and fill in the real invocation
 * (both take PLINK-format genotype data, similar to the king_plink
 * module's conversion step).
 */
process READ_KIN {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.read_kin.stub"), emit: result

    script:
    """
    echo "READ/KIN are not wired up yet: ${meta.reason}" >&2
    exit 1
    """
}
