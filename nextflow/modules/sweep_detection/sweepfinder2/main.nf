/*
 * SweepFinder2
 * Status: runnable (bioconda: sweepfinder2, cyvcf2)
 *
 * Converts the VCF to SweepFinder2's .freq format (bin/vcf_to_sweepfinder_freq.py,
 * unfolded - see that script's ancestral-allele assumption), builds the
 * background/empirical spectrum with `-f`, then scans with `-l G` across
 * the whole contig. `args` must include both the scan mode and its grid
 * argument; challenge configs default to `-l 1000` (1,000 scan positions).
 */
process SWEEPFINDER2 {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/sweepfinder2" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.sf2.out"), emit: result

    script:
    """
    vcf_to_sweepfinder_freq.py ${input_file} ${meta.challenge}.freq
    SweepFinder2 -f ${meta.challenge}.freq ${meta.challenge}.spect
    SweepFinder2 ${meta.args} ${meta.challenge}.freq ${meta.challenge}.spect ${meta.challenge}.sf2.out
    """
}
