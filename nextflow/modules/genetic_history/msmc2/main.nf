/*
 * MSMC2 (Multiple Sequentially Markovian Coalescent)
 * Status: runnable (bioconda: msmc2, cyvcf2) for 'vcf' input;
 * runnable-unverified (bioconda: msmc2, bcftools, samtools, bedtools,
 * cyvcf2) for 'bam' input - see module header below, untested pending
 * the real BAM archives (needs Synapse challenge-team auth to download -
 * see conf/genetic_history/bottleneck_low_coverage.config).
 *
 * Two entry points depending on the challenge's `input_type`:
 *   - vcf: builds a multihetsep file directly from genotypes
 *          (bin/vcf_to_multihetsep.py, no mask - every site treated as
 *          callable, since no coverage/reference info exists for these
 *          challenges at all).
 *   - bam: the GHIST *_low_coverage challenges' real input - an archive
 *          of per-individual BAMs + popfile.txt + reference FASTA.
 *          bin/bam_archive_to_merged_vcf.sh picks up to 2 individuals
 *          (one per population if `split_isolation_low_coverage`'s
 *          popfile lists 2 populations - giving MSMC2 the cross-
 *          population haplotype pairing its split-time estimate needs;
 *          otherwise the first 2 individuals from the single population,
 *          for `bottleneck_low_coverage`), calls each with bcftools
 *          mpileup/call, builds a real per-sample depth-based
 *          callability mask (bin/bam_to_depth_mask.py), and intersects
 *          them - then bin/vcf_to_multihetsep.py --mask uses that real
 *          mask instead of the "every site is callable" assumption the
 *          vcf path has to make.
 */
process MSMC2 {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/msmc2" }, mode: 'copy'
    conda { meta.conda }

    input:
    // ref_file is a harmless placeholder for 'vcf' challenges (that
    // script branch never reads it) - stageAs forces a distinct staged
    // name so it doesn't collide with main_file when both resolve to
    // the same underlying path (see workflows/genetic_history.nf).
    tuple val(meta), path(main_file), path(ref_file, stageAs: 'ref_placeholder.*')

    output:
    tuple val(meta), path("${meta.challenge}.msmc2.final.txt"), emit: result

    script:
    // bioconda's msmc2 recipe installs the binary as 'msmc2_Linux' (or
    // 'msmc2_Darwin' on macOS) - there is no plain 'msmc2' alias -
    // confirmed when this module was first run against real data, so
    // resolve the actual binary name at runtime instead of hardcoding it.
    def msmc2Bin = 'MSMC2_BIN=$(compgen -c | grep -E "^msmc2(_Linux|_Darwin)?$" | head -1)'
    if (meta.input_type == 'vcf')
        """
        ${msmc2Bin}
        vcf_to_multihetsep.py ${main_file} ${meta.challenge}.multihetsep.txt
        "\$MSMC2_BIN" ${meta.args} -o ${meta.challenge}.msmc2 ${meta.challenge}.multihetsep.txt
        """
    else
        """
        ${msmc2Bin}
        bam_archive_to_merged_vcf.sh ${main_file} ${ref_file} ${meta.challenge}
        vcf_to_multihetsep.py ${meta.challenge}.merged.vcf.gz ${meta.challenge}.multihetsep.txt \
            --mask ${meta.challenge}.mask.bed --max-samples 2
        "\$MSMC2_BIN" ${meta.args} -o ${meta.challenge}.msmc2 ${meta.challenge}.multihetsep.txt
        """
}
