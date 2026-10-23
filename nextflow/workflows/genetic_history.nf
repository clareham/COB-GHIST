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
 *
 * 'vcf' challenges first pass through NORMALISE_VCF (biallelic SNPs only,
 * repaired header, excluded-sites BED); 'bam' challenges (the
 * *_low_coverage ones) ship an archive of per-individual BAMs plus the
 * reference FASTA and are only supported by the PSMC/MSMC2 BAM paths.
 */

include { NORMALISE_VCF    } from '../modules/common/normalise_vcf/main.nf'
include { PSMC             } from '../modules/genetic_history/psmc/main.nf'
include { MSMC2            } from '../modules/genetic_history/msmc2/main.nf'
include { FASTSIMCOAL2     } from '../modules/genetic_history/fastsimcoal2/main.nf'
include { MRPAST           } from '../modules/genetic_history/mrpast/main.nf'
include { STAIRWAY_PLOT2   } from '../modules/genetic_history/stairway_plot2/main.nf'
include { PHLASH           } from '../modules/genetic_history/phlash/main.nf'
include { FASTGLOBETROTTER } from '../modules/genetic_history/fastglobetrotter/main.nf'
include { ARCHIE2          } from '../modules/genetic_history/archie2/main.nf'

def toolMeta(String chName, Map chConf, String _toolName, Map toolConf) {
    return [
        category  : 'genetic_history',
        challenge : chName,
        input_type: chConf.input_type,
        conda     : toolConf.conda,
        args      : toolConf.args ?: '',
        options   : toolConf.options ?: [:],
        reason    : toolConf.reason ?: '',
    ]
}

// Entries of the (tool, meta, vcf, tbi, excluded_bed) fan-out channel for one tool.
def forTool(byTool, String name) {
    return byTool.filter { entry -> entry[0] == name }
}

workflow GENETIC_HISTORY {
    main:
    def challenges = params.genetic_history.challenges

    def vcfInputs = []
    def bamByTool = [:].withDefault { [] }

    challenges.each { chName, chConf ->
        if (!chConf.enabled) return
        def enabledTools = chConf.tools.findAll { _name, toolConf -> toolConf.enabled }
        if (!enabledTools) return

        def mainName = chConf.input_type == 'vcf' ? chConf.testing_vcf_name : chConf.testing_archive_name
        def mainFile = GhistUtils.testingFile(params.data_dir, mainName)
        if (!mainFile.exists()) {
            log.warn "genetic_history/${chName}: testing file not found (${mainFile}), skipping challenge"
            return
        }

        if (chConf.input_type == 'vcf') {
            vcfInputs << tuple([category: 'genetic_history', challenge: chName], file(mainFile))
            return
        }

        def refFile = GhistUtils.testingFile(params.data_dir, chConf.testing_fasta_name)
        if (!refFile.exists()) {
            log.warn "genetic_history/${chName}: reference FASTA not found (${refFile}), skipping challenge"
            return
        }
        enabledTools.each { toolName, toolConf ->
            if (!(toolName in ['psmc', 'msmc2'])) {
                log.warn "genetic_history/${chName}: ${toolName} has no BAM workflow, skipping"
                return
            }
            // The empty list leaves the excluded-sites BED input unstaged.
            bamByTool[toolName] << tuple(toolMeta(chName, chConf, toolName, toolConf), file(mainFile), file(refFile), [])
        }
    }

    NORMALISE_VCF(channel.fromList(vcfInputs))

    // One (tool, meta, vcf, tbi, excluded_bed) entry per enabled tool of
    // each normalised challenge.
    def byTool = NORMALISE_VCF.out.vcf.flatMap { chMeta, vcf, tbi, excluded ->
        def chConf = challenges[chMeta.challenge]
        chConf.tools.findAll { _name, toolConf -> toolConf.enabled }.collect { toolName, toolConf ->
            tuple(toolName, toolMeta(chMeta.challenge, chConf, toolName, toolConf), vcf, tbi, excluded)
        }
    }
    PSMC(
        forTool(byTool, 'psmc').map { _tool, meta, vcf, _tbi, excluded -> tuple(meta, vcf, vcf, excluded) }
            .mix(channel.fromList(bamByTool['psmc']))
    )
    MSMC2(
        forTool(byTool, 'msmc2').map { _tool, meta, vcf, _tbi, excluded -> tuple(meta, vcf, vcf, excluded) }
            .mix(channel.fromList(bamByTool['msmc2']))
    )
    FASTSIMCOAL2(forTool(byTool, 'fastsimcoal2').map { _tool, meta, vcf, _tbi, excluded -> tuple(meta, vcf, excluded) })
    MRPAST(forTool(byTool, 'mrpast').map { _tool, meta, vcf, _tbi, _excluded -> tuple(meta, vcf) })
    STAIRWAY_PLOT2(forTool(byTool, 'stairway_plot2').map { _tool, meta, vcf, _tbi, _excluded -> tuple(meta, vcf) })
    PHLASH(forTool(byTool, 'phlash').map { _tool, meta, vcf, _tbi, _excluded -> tuple(meta, vcf) })
    FASTGLOBETROTTER(forTool(byTool, 'fastglobetrotter').map { _tool, meta, vcf, _tbi, _excluded -> tuple(meta, vcf) })
    ARCHIE2(forTool(byTool, 'archie2').map { _tool, meta, vcf, _tbi, _excluded -> tuple(meta, vcf) })

    emit:
    psmc         = PSMC.out.result
    msmc2        = MSMC2.out.result
    fastsimcoal2 = FASTSIMCOAL2.out.result
    mrpast       = MRPAST.out.result
}
