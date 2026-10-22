/*
 * Genetic History category workflow.
 *
 * Fans out over `params.genetic_history.challenges` (level 2 of the
 * config hierarchy) and, within each enabled challenge, over its
 * `tools {}` map (level 3), routing each enabled (challenge, tool) pair
 * into that tool's shared process. A tool with no enabled challenges in
 * this run simply gets an empty input channel and never executes - this
 * is how the many `status: 'stub'`/`'excluded'` entries stay inert by
 * default without extra conditional logic.
 */

include { PSMC             } from '../modules/genetic_history/psmc/main.nf'
include { MSMC2            } from '../modules/genetic_history/msmc2/main.nf'
include { FASTSIMCOAL2     } from '../modules/genetic_history/fastsimcoal2/main.nf'
include { MRPAST           } from '../modules/genetic_history/mrpast/main.nf'
include { STAIRWAY_PLOT2   } from '../modules/genetic_history/stairway_plot2/main.nf'
include { PHLASH           } from '../modules/genetic_history/phlash/main.nf'
include { FASTGLOBETROTTER } from '../modules/genetic_history/fastglobetrotter/main.nf'
include { ARCHIE2          } from '../modules/genetic_history/archie2/main.nf'

workflow GENETIC_HISTORY {
    main:
    def byTool = [:].withDefault { [] }

    params.genetic_history.challenges.each { chName, chConf ->
        if (!chConf.enabled) return

        // 'vcf' challenges: one file, no separate reference exists.
        // 'bam' challenges (the *_low_coverage ones): an archive of
        // per-individual BAMs + popfile.txt, plus the reference FASTA
        // needed to call variants from them.
        def mainName = chConf.input_type == 'vcf' ? chConf.testing_vcf_name : chConf.testing_archive_name
        def mainFile = GhistUtils.testingFile(params.data_dir, mainName)
        if (!mainFile.exists()) {
            log.warn "genetic_history/${chName}: testing file not found (${mainFile}), skipping challenge"
            return
        }
        // PSMC/MSMC2 take (main_file, ref_file) uniformly across both
        // input types - 'vcf' challenges pass a harmless placeholder
        // (reusing mainFile) since that script branch never reads it.
        def refFile = mainFile
        if (chConf.input_type == 'bam') {
            refFile = GhistUtils.testingFile(params.data_dir, chConf.testing_fasta_name)
            if (!refFile.exists()) {
                log.warn "genetic_history/${chName}: reference FASTA not found (${refFile}), skipping challenge"
                return
            }
        }

        chConf.tools.each { toolName, toolConf ->
            if (!toolConf.enabled) return
            def meta = [
                category  : 'genetic_history',
                challenge : chName,
                input_type: chConf.input_type,
                conda     : toolConf.conda,
                args      : toolConf.args ?: '',
                options   : toolConf.options ?: [:],
                reason    : toolConf.reason ?: '',
            ]
            if (toolName in ['psmc', 'msmc2']) {
                byTool[toolName] << tuple(meta, file(mainFile), file(refFile))
            } else {
                byTool[toolName] << tuple(meta, file(mainFile))
            }
        }
    }

    PSMC(channel.fromList(byTool['psmc']))
    MSMC2(channel.fromList(byTool['msmc2']))
    FASTSIMCOAL2(channel.fromList(byTool['fastsimcoal2']))
    MRPAST(channel.fromList(byTool['mrpast']))
    STAIRWAY_PLOT2(channel.fromList(byTool['stairway_plot2']))
    PHLASH(channel.fromList(byTool['phlash']))
    FASTGLOBETROTTER(channel.fromList(byTool['fastglobetrotter']))
    ARCHIE2(channel.fromList(byTool['archie2']))

    emit:
    psmc         = PSMC.out.result
    msmc2        = MSMC2.out.result
    fastsimcoal2 = FASTSIMCOAL2.out.result
    mrpast       = MRPAST.out.result
}
