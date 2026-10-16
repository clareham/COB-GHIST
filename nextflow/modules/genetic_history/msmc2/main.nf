/*
 * MSMC2 survey: one input per contig, explicit phase ambiguity, and
 * separate within/cross-population runs with relative cross-coalescence.
 * Mutation-scaled rates are preserved; time in rccr.tsv uses the configured
 * mutation rate. These curves are not a fitted challenge submission.
 */
process MSMC2 {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/msmc2" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(main_file), path(ref_file, stageAs: 'ref_placeholder.*')

    output:
    tuple val(meta), path("${meta.challenge}.msmc2.final.txt"), emit: result
    tuple val(meta), path("${meta.challenge}.msmc2.run"), emit: details

    script:
    def options = meta.options + [threads: task.cpus, extra_args: meta.args]
    def optionsJson = groovy.json.JsonOutput.toJson(options)
    def prepare = meta.input_type == 'vcf' ? '' :
        "bam_archive_to_merged_vcf.sh ${main_file} ${ref_file} ${meta.challenge} ${options.min_depth} ${options.max_depth}"
    def inputVcf = meta.input_type == 'vcf' ? main_file : "${meta.challenge}.merged.vcf.gz"
    def bamArgs = meta.input_type == 'vcf' ? '' :
        "--mask ${meta.challenge}.mask.bed --population-map ${meta.challenge}.populations.tsv"
    """
    cat > msmc_options.json <<'JSON'
    ${optionsJson}
    JSON
    ${prepare}
    msmc_survey.py ${inputVcf} ${meta.challenge} msmc_options.json ${bamArgs}
    """
}
