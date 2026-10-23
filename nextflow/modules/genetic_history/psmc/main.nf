/*
 * PSMC survey for one diploid sample. Simulated variant-only VCFs assume
 * declared contigs are callable; missing/invalid genotypes mask windows.
 * BAMs use configured quality/depth filters, not a reference consensus
 * masquerading as a diploid sequence. PSMC does not fit population splits.
 */
process PSMC {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/psmc" }, mode: 'copy'
    conda { meta.conda }

    input:
    // VCF challenges: (normalised VCF, VCF again as a harmless ref placeholder,
    // excluded-sites BED). BAM challenges: (archive, reference FASTA, []).
    tuple val(meta), path(main_file), path(ref_file, stageAs: 'ref_placeholder.*'), path(excluded_bed)

    output:
    tuple val(meta), path("${meta.challenge}.psmc"), emit: result
    tuple val(meta), path("${meta.challenge}.psmc.assumptions.json"), emit: assumptions

    script:
    def optionsJson = groovy.json.JsonOutput.toJson(meta.options + [input_type: meta.input_type, args: meta.args])
    def prepare
    if (meta.input_type == 'vcf')
        prepare = """
        vcf_sample_to_psmcfa.py ${main_file} ${meta.challenge}.psmcfa ${meta.options.assume_callable ? '--assume-callable' : ''} ${excluded_bed ? "--exclude ${excluded_bed}" : ''}
        """
    else
        prepare = """
        bam_archive_samples.py ${main_file} extracted ${meta.challenge}.populations.tsv --max-samples 1
        BAM=\$(tail -n +2 ${meta.challenge}.populations.tsv | cut -f3)
        bam_sample_to_psmcfa.sh ${ref_file} "\$BAM" ${meta.challenge} ${meta.options.min_depth} ${meta.options.max_depth}
        """
    """
    cat > ${meta.challenge}.psmc.assumptions.json <<'JSON'
    ${optionsJson}
    JSON
    ${prepare}
    if ! awk '!/^>/ && /K/ {found=1} END {exit !found}' ${meta.challenge}.psmcfa; then
        echo "ERROR: PSMC input has no callable heterozygous windows" >&2
        exit 1
    fi
    psmc ${meta.args} -o ${meta.challenge}.psmc ${meta.challenge}.psmcfa
    if [ ! -s ${meta.challenge}.psmc ]; then
        echo "ERROR: PSMC did not produce a nonempty result" >&2
        exit 1
    fi
    """
}
