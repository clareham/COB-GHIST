/*
 * SHAPEIT5
 * Status: excluded by default (testing VCFs are already phased).
 *
 * Optional rephasing only, not a relatedness estimator. Bioconda installs
 * SHAPEIT5_phase_common. Each populated contig is phased separately and
 * concatenated; without a genetic map SHAPEIT5 assumes 1 cM/Mb.
 */
process SHAPEIT5 {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/shapeit5" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.phased.bcf"), emit: result

    script:
    """
    bgzip -c ${input_file} > ${meta.challenge}.vcf.gz
    tabix -p vcf ${meta.challenge}.vcf.gz
    bcftools index -s ${meta.challenge}.vcf.gz > contig_stats.tsv
    awk '\$3 > 0 {print \$1}' contig_stats.tsv > contigs.txt
    if [ ! -s contigs.txt ]; then
        echo "ERROR: no populated contigs found for SHAPEIT5" >&2
        exit 1
    fi
    mkdir phased
    : > phased.list
    i=0
    while IFS= read -r region; do
        SHAPEIT5_phase_common ${meta.args} --thread ${task.cpus} --input ${meta.challenge}.vcf.gz --region "\$region" --output "phased/\$i.bcf"
        bcftools index "phased/\$i.bcf"
        printf '%s\\n' "phased/\$i.bcf" >> phased.list
        i=\$((i + 1))
    done < contigs.txt
    bcftools concat -f phased.list -Ob -o ${meta.challenge}.phased.bcf
    if [ ! -s ${meta.challenge}.phased.bcf ]; then
        echo "ERROR: SHAPEIT5 did not produce a nonempty phased BCF" >&2
        exit 1
    fi
    """
}
