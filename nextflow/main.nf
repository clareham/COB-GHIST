#!/usr/bin/env nextflow
/*
 * GHIST 2026 survey pipeline - entry point.
 *
 * See nextflow.config for the config hierarchy (category -> challenge ->
 * tool) and docs/existingtools.md (in the parent repo) for the tool
 * catalogue this scaffold is built from.
 *
 * Only *.testing.* / *_testing.* files are ever read - see
 * lib/GhistUtils.groovy. *.final.* sets are the held-out GHIST scoring
 * data and must never be consumed by this pipeline.
 */

nextflow.enable.dsl = 2

include { GENETIC_HISTORY } from './workflows/genetic_history.nf'
include { SWEEP_DETECTION } from './workflows/sweep_detection.nf'
include { RELATEDNESS     } from './workflows/relatedness.nf'

workflow {
    log.info """
    GHIST 2026 survey pipeline
    ---------------------------
    data_dir : ${params.data_dir}
    outdir   : ${params.outdir}
    genetic_history : ${params.run_genetic_history}
    sweep_detection : ${params.run_sweep_detection}
    relatedness     : ${params.run_relatedness}
    """.stripIndent()

    if (params.run_genetic_history) GENETIC_HISTORY()
    if (params.run_sweep_detection) SWEEP_DETECTION()
    if (params.run_relatedness)     RELATEDNESS()
}
