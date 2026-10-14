/*
 * Stairway Plot 2
 * Status: stub - distributed as a Java .jar via GitHub release, no
 * conda/PyPI package. To activate: download the release jar, set
 * `params.<category>.challenges.<challenge>.tools.stairway_plot2.jar`
 * to its path, flip `enabled: true` in the relevant challenge config,
 * and fill in the real invocation below (it needs a 2-pop/1-pop blueprint
 * file built from bin/vcf_to_sfs.py's SFS, which Stairway Plot 2 expects
 * in its own blueprint text format, not fastsimcoal2's .obs format).
 */
process STAIRWAY_PLOT2 {
    tag "${meta.challenge}"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.stairway.stub"), emit: result

    script:
    """
    echo "Stairway Plot 2 is not wired up yet: ${meta.reason}" >&2
    exit 1
    """
}
