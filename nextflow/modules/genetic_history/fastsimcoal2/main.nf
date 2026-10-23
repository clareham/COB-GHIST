process FASTSIMCOAL2 {
    tag "${meta.challenge}"

    conda "${meta.conda}"
    cpus { (meta.options?.threads ?: 1) as int }
    // 100k simulations x 40 ECM cycles on a 100+ haplotype SFS, plus the
    // pure-Python SFS build over up to 5M records, exceeds the 2 h default.
    time '8.h'

    input:
    tuple val(meta), path(input_file), path(excluded_bed)

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
        --options '${optionsJson}' \\
        --exclude-bed "${excluded_bed}"
    cd "\$outdir"

    # fsc27093 is the executable shipped by Bioconda's fastsimcoal2=27093.
    # -m reads the folded (minor-allele) *_MAFpop0.obs spectrum and -d the
    # unfolded *_jointDAFpop1_0.obs spectrum; --foldedSFS is a simulation
    # output option and does not select the observed SFS type.
    if [[ "${meta.challenge}" == "bottleneck" || "${meta.challenge}" == "growth" ]]; then
        fsc27093 -t "${meta.challenge}.tpl" -e "${meta.challenge}.est" \\
            -n ${simulations} -M -L ${cycles} -m -q -c ${task.cpus}
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
