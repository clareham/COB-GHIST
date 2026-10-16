process MRPAST {
    tag "${meta.challenge}"

    conda "${meta.conda}"
    cpus { (meta.options?.jobs ?: 1) as int }

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
