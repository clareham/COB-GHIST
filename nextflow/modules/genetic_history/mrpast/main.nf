process MRPAST {
    tag "${meta.challenge}"

    conda "${meta.conda}"
    cpus { (meta.options?.jobs ?: 1) as int }
    // tsinfer + tsdate over a 100 Mb contig with >100 haplotypes takes
    // hours; the 2 h default would kill the task mid-inference.
    time '8.h'
    memory '16.GB'

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.mrpast_out"), emit: result

    script:
    def optionsJson = groovy.json.JsonOutput.toJson(meta.options ?: [:])
    """
    set -euo pipefail
    outdir="${meta.challenge}.mrpast_out"
    mkdir -p "\$outdir"

    mrpast_survey.py \\
        --vcf "${input_file}" \\
        --challenge "${meta.challenge}" \\
        --outdir "\$outdir" \\
        --options '${optionsJson}'
    """
}
