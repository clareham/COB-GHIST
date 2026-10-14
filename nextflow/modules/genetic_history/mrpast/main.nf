/*
 * mrpast (demographic inference from Ancestral Recombination Graphs)
 * Status: runnable-unverified (bioconda: mrpast)
 *
 * CAVEAT: mrpast is ARG-based - its natural input is a tree sequence
 * (.trees), not a raw VCF. This module takes the pragmatic survey-stage
 * shortcut of inferring an ARG from the VCF with tsinfer first (also
 * bioconda-packaged), rather than requiring a separate ARG-building
 * pipeline up front. We have not independently verified mrpast's exact
 * CLI flags beyond `--help` - treat this as a starting point to check
 * against the real mrpast docs (https://aprilweilab.github.io/) before
 * trusting results.
 */
process MRPAST {
    tag "${meta.challenge}"
    publishDir { "${params.outdir}/${meta.category}/${meta.challenge}/mrpast" }, mode: 'copy'
    conda { meta.conda }

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("${meta.challenge}.mrpast_out"), emit: result

    script:
    """
    bgzip -c ${input_file} > ${meta.challenge}.vcf.gz
    tabix -p vcf ${meta.challenge}.vcf.gz

    python3 - <<PY
import tsinfer, cyvcf2, numpy as np

vcf = cyvcf2.VCF("${meta.challenge}.vcf.gz")
with tsinfer.SampleData(path="${meta.challenge}.samples", sequence_length=1e9) as samples:
    for s in vcf.samples:
        samples.add_individual(ploidy=2, metadata={"name": s})
    for rec in vcf:
        samples.add_site(rec.POS, [g for gt in rec.genotypes for g in gt[:2]])

ts = tsinfer.infer(tsinfer.load("${meta.challenge}.samples"))
ts.dump("${meta.challenge}.trees")
PY

    mkdir -p ${meta.challenge}.mrpast_out
    mrpast ${meta.args} --trees ${meta.challenge}.trees --out-dir ${meta.challenge}.mrpast_out \
        || echo "mrpast invocation failed/unverified - check CLI flags against 'mrpast --help'" >&2
    """
}
