/*
 * FlexSweep 2.0 survey: simulate -> features -> train -> VCF features -> predict.
 * No pretrained model is assumed. Training is expensive; survey-sized training
 * is exploratory, not a calibrated demographic/background-selection model.
 * CLI contracts and version-specific repairs are documented in the helper.
 */
process FLEXSWEEP {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/flexsweep" }, mode: 'copy'
    conda { meta.conda }
    cpus { (meta.options?.threads ?: 4) as int }
    memory '16 GB'
    time '24h'

    input:
    tuple val(meta), path(input_file), path(accessible_bed)

    output:
    tuple val(meta), path("${meta.challenge}.flexsweep"), emit: result

    script:
    if (meta.args?.trim()) {
        error 'FlexSweep survey uses structured tool options, not unvalidated args.'
    }
    def optionsJson = groovy.json.JsonOutput.toJson(meta.options ?: [:])
    def quotedInput = "'" + input_file.toString().replace("'", "'\\''") + "'"
    def quotedMask = "'" + accessible_bed.toString().replace("'", "'\\''") + "'"
    """
    cat > flexsweep-options.json <<'FLEXSWEEP_OPTIONS'
    ${optionsJson}
    FLEXSWEEP_OPTIONS
    flexsweep_survey.py --vcf ${quotedInput} --mask ${quotedMask} \\
        --options flexsweep-options.json --threads ${task.cpus} \\
        --output '${meta.challenge}.flexsweep'
    """
}
