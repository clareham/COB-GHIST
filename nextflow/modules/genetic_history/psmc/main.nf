/*
 * PSMC (Pairwise Sequentially Markovian Coalescent)
 * Status: runnable (bioconda: psmc, bcftools, cyvcf2) for 'vcf' input;
 * runnable-unverified (bioconda: psmc, bcftools, samtools) for 'bam'
 * input - see module header below, untested pending the real BAM
 * archives (needs Synapse challenge-team auth to download - see
 * conf/genetic_history/bottleneck_low_coverage.config).
 *
 * Two entry points depending on the challenge's `input_type`:
 *   - vcf: derive a PSMC-fa call string directly from one sample's
 *          genotypes (bin/vcf_sample_to_psmcfa.py). No reference FASTA
 *          is needed/available for these simulated contigs, so the
 *          classic `bcftools consensus | vcfutils.pl vcf2fq` route is
 *          skipped in favour of building the T/K/N string straight
 *          from GT calls.
 *   - bam: the GHIST *_low_coverage challenges hand out an archive of
 *          per-individual BAMs (~2-4x coverage) + a `popfile.txt`
 *          (sample -> population) + a reference FASTA - this IS PSMC's
 *          real intended input. Extracts the archive, picks the first
 *          individual listed in popfile.txt (consistent with how other
 *          multi-sample challenges pick one representative sample for
 *          this single-individual method), and runs the classic
 *          bcftools mpileup | bcftools call | vcfutils.pl vcf2fq |
 *          fq2psmcfa pipeline (bin/bam_sample_to_psmcfa.sh) - PSMC's own
 *          documented usage. For `split_isolation_low_coverage` (2
 *          populations), PSMC only ever reports a single-population
 *          curve for whichever individual it picks; see MSMC2 for the
 *          cross-population comparison this challenge actually wants.
 */
process PSMC {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/psmc" }, mode: 'copy'
    conda { meta.conda }

    input:
    // ref_file is a harmless placeholder for 'vcf' challenges (that
    // script branch never reads it) - stageAs forces a distinct staged
    // name so it doesn't collide with main_file when both resolve to
    // the same underlying path (see workflows/genetic_history.nf).
    tuple val(meta), path(main_file), path(ref_file, stageAs: 'ref_placeholder.*')

    output:
    tuple val(meta), path("${meta.challenge}.psmc"), emit: result

    script:
    if (meta.input_type == 'vcf')
        """
        vcf_sample_to_psmcfa.py ${main_file} ${meta.challenge}.psmcfa
        psmc ${meta.args} -o ${meta.challenge}.psmc ${meta.challenge}.psmcfa
        """
    else
        """
        mkdir -p extracted
        tar xzf ${main_file} -C extracted
        ARCHIVE_DIR=\$(find extracted -mindepth 1 -maxdepth 1 -type d | head -1)
        POPFILE="\${ARCHIVE_DIR}/popfile.txt"
        # popfile.txt has a header row ('sample\\tpop') - skip it to get the
        # first real sample, not the literal header text.
        SAMPLE=\$(tail -n +2 "\$POPFILE" | head -1 | cut -f1)

        BAM=""
        for f in "\${ARCHIVE_DIR}"/*.bam; do
            SM=\$(samtools view -H "\$f" | awk -F'\\t' '/^@RG/{for(i=1;i<=NF;i++) if(\$i ~ /^SM:/) print substr(\$i,4)}' | head -1)
            if [ "\$SM" = "\$SAMPLE" ]; then BAM="\$f"; break; fi
        done
        if [ -z "\$BAM" ]; then
            echo "ERROR: no BAM in ${main_file} has read-group SM matching popfile sample '\$SAMPLE'" >&2
            exit 1
        fi

        bam_sample_to_psmcfa.sh ${ref_file} "\$BAM" ${meta.challenge}
        psmc ${meta.args} -o ${meta.challenge}.psmc ${meta.challenge}.psmcfa
        """
}
