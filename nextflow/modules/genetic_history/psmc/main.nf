/*
 * PSMC (Pairwise Sequentially Markovian Coalescent)
 * Status: runnable (bioconda: psmc, cyvcf2)
 *
 * Two entry points depending on the challenge's `input_type`:
 *   - vcf:   derive a PSMC-fa call string directly from one sample's
 *            genotypes (bin/vcf_sample_to_psmcfa.py). No reference FASTA
 *            is needed/available for these simulated contigs, so the
 *            classic `bcftools consensus | vcfutils.pl vcf2fq` route is
 *            skipped in favour of building the T/K/N string straight
 *            from GT calls.
 *   - fasta: the GHIST *_low_coverage challenges hand out a single
 *            haploid consensus sequence (no second haplotype, no
 *            heterozygosity signal). CAVEAT: classic PSMC needs per-site
 *            het/hom calls from a single *diploid* genome; a lone haploid
 *            sequence has none, so bin/fasta_to_pseudo_psmcfa.py produces
 *            a flat, uninformative curve. Treat that path's output as a
 *            structural placeholder, not a validated inference.
 */
process PSMC {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/psmc" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.psmc"), emit: result

    script:
    // bin/ scripts are on PATH automatically (Nextflow convention: <project>/bin/).
    def converter = meta.input_type == 'vcf' ?
        "vcf_sample_to_psmcfa.py ${input_file} ${meta.challenge}.psmcfa" :
        "fasta_to_pseudo_psmcfa.py ${input_file} ${meta.challenge}.psmcfa"
    """
    ${converter}
    psmc ${meta.args} -o ${meta.challenge}.psmc ${meta.challenge}.psmcfa
    """
}
