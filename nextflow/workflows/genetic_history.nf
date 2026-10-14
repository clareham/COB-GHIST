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
        def fileName = chConf.input_type == 'vcf' ? chConf.testing_vcf_name : chConf.testing_fasta_name
        def inputFile = GhistUtils.testingFile(params.data_dir, fileName)
        if (!inputFile.exists()) {
            log.warn "genetic_history/${chName}: testing file not found (${inputFile}), skipping challenge"
            return
        }
        chConf.tools.each { toolName, toolConf ->
            if (!toolConf.enabled) return
            def meta = [
                category  : 'genetic_history',
                challenge : chName,
                input_type: chConf.input_type,
                conda     : toolConf.conda,
                args      : toolConf.args ?: '',
                reason    : toolConf.reason ?: '',
            ]
            byTool[toolName] << tuple(meta, file(inputFile))
        }
    }

    PSMC(Channel.fromList(byTool['psmc']))
    MSMC2(Channel.fromList(byTool['msmc2']))
    FASTSIMCOAL2(Channel.fromList(byTool['fastsimcoal2']))
    MRPAST(Channel.fromList(byTool['mrpast']))
    STAIRWAY_PLOT2(Channel.fromList(byTool['stairway_plot2']))
    PHLASH(Channel.fromList(byTool['phlash']))
    FASTGLOBETROTTER(Channel.fromList(byTool['fastglobetrotter']))
    ARCHIE2(Channel.fromList(byTool['archie2']))

    emit:
    psmc         = PSMC.out.result
    msmc2        = MSMC2.out.result
    fastsimcoal2 = FASTSIMCOAL2.out.result
    mrpast       = MRPAST.out.result
}
