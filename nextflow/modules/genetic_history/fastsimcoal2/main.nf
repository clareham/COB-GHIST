process FASTSIMCOAL2 {
    tag "${meta.challenge}"

    conda "${meta.conda}"
    cpus { (meta.options?.threads ?: 1) as int }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.fastsimcoal2_out"), emit: result

    script:
    def optionsJson = groovy.json.JsonOutput.toJson(meta.options ?: [:])
    def simulations = (meta.options?.simulations ?: 100000) as int
    def cycles = (meta.options?.ecm_cycles ?: 40) as int
    """
    set -euo pipefail
    outdir="${meta.challenge}.fastsimcoal2_out"
    mkdir -p "\$outdir"
    fastsimcoal2_survey.py \\
        --vcf "${input_file}" \\
        --challenge "${meta.challenge}" \\
        --outdir "\$outdir" \\
        --options '${optionsJson}'
    cd "\$outdir"

    # fsc27093 is the executable shipped by Bioconda's fastsimcoal2=27093.
    if [[ "${meta.challenge}" == "bottleneck" || "${meta.challenge}" == "growth" ]]; then
        fsc27093 -t "${meta.challenge}.tpl" -e "${meta.challenge}.est" \\
            -n ${simulations} -M -L ${cycles} --foldedSFS -q -c ${task.cpus}
    else
        fsc27093 -t "${meta.challenge}.tpl" -e "${meta.challenge}.est" \\
            -n ${simulations} -M -L ${cycles} -d -q -c ${task.cpus}
    fi
    result_dir="${meta.challenge}"
    if [[ ! -d "\$result_dir" ]] || ! find "\$result_dir" -type f -name '*.bestlhoods' -print -quit | grep -q .; then
        echo "fastsimcoal2 completed without producing a best-likelihood result" >&2
        exit 1
    fi
    """
}
