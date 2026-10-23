/*
 * SweepFinder2
 * Status: runnable (bioconda: sweepfinder2, cyvcf2)
 *
 * Converts the VCF to SweepFinder2's .freq format (bin/vcf_to_sweepfinder_freq.py,
 * folded by default because these VCFs lack ancestral-allele tags), builds the
 * background/empirical spectrum with `-f`, then scans with `-l G` across
 * the whole contig. `args` must include both the scan mode and its grid
 * argument; challenge configs default to `-l 1000` (1,000 scan positions).
 */
process SWEEPFINDER2 {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/sweepfinder2" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file), path(accessible_bed)

    output:
    tuple val(meta), path("${meta.challenge}.sf2.out"), emit: result

    script:
    """
    vcf_to_sweepfinder_freq.py ${input_file} ${meta.challenge}.freq --mask ${accessible_bed} ${meta.options.folded ? '--folded' : ''}
    SweepFinder2 -f ${meta.challenge}.freq ${meta.challenge}.spect
    SweepFinder2 ${meta.args} ${meta.challenge}.freq ${meta.challenge}.spect ${meta.challenge}.sf2.out
    if [ ! -s ${meta.challenge}.sf2.out ]; then
        echo "ERROR: SweepFinder2 did not produce a nonempty scan" >&2
        exit 1
    fi
    """
}
