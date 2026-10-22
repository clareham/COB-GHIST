/*
 * Relatedness Estimation category workflow. Same fan-out pattern as
 * workflows/genetic_history.nf - see that file's header comment. VCF
 * challenges pass through NORMALISE_VCF first; this also guarantees the
 * biallelic-only input PLINK2's --make-bed requires.
 *
 * relatedness1 contributes no entries here: its challenge config has an
 * empty `tools {}` map (data/tool format mismatch - see
 * conf/relatedness/relatedness1.config).
 */

include { NORMALISE_VCF } from '../modules/common/normalise_vcf/main.nf'
include { KING_PLINK    } from '../modules/relatedness/king_plink/main.nf'
include { NGSRELATE     } from '../modules/relatedness/ngsrelate/main.nf'
include { ADMIXTURE     } from '../modules/relatedness/admixture/main.nf'
include { SHAPEIT5      } from '../modules/relatedness/shapeit5/main.nf'
include { READ_KIN      } from '../modules/relatedness/read_kin/main.nf'
include { ICMLKIN       } from '../modules/relatedness/icmlkin/main.nf'
include { GRUPS_RS      } from '../modules/relatedness/grups_rs/main.nf'

def toolMeta(String chName, Map chConf, String _toolName, Map toolConf) {
    return [
        category  : 'relatedness',
        challenge : chName,
        input_type: chConf.input_type,
        autosomes : chConf.autosomes,
        conda     : toolConf.conda,
        args      : toolConf.args ?: '',
        k         : toolConf.k ?: 3,
        reason    : toolConf.reason ?: '',
    ]
}

// (meta, vcf) entries of the (tool, meta, vcf, tbi, excluded_bed) fan-out channel for one tool.
def forTool(byTool, String name) {
    return byTool.filter { entry -> entry[0] == name }
        .map { _tool, meta, vcf, _tbi, _excluded -> tuple(meta, vcf) }
}

workflow RELATEDNESS {
    main:
    def challenges = params.relatedness.challenges

    def vcfInputs = []
    challenges.each { chName, chConf ->
        if (!chConf.enabled) return
        def enabledTools = chConf.tools.findAll { _name, toolConf -> toolConf.enabled }
        if (!enabledTools) return
        def fileName = chConf.input_type == 'vcf' ? chConf.testing_vcf_name : chConf.testing_csv_name
        def inputFile = GhistUtils.testingFile(params.data_dir, fileName)
        if (!inputFile.exists()) {
            log.warn "relatedness/${chName}: testing file not found (${inputFile}), skipping challenge"
            return
        }
        if (chConf.input_type != 'vcf') {
            log.warn "relatedness/${chName}: no tool accepts ${chConf.input_type} input, skipping challenge"
            return
        }
        enabledTools.each { toolName, _toolConf ->
            if (toolName in ['king_plink', 'admixture'] &&
                (!(chConf.autosomes instanceof Integer) || chConf.autosomes < 1)) {
                error "relatedness/${chName}: ${toolName} requires a positive diploid autosome count"
            }
        }
        vcfInputs << tuple([category: 'relatedness', challenge: chName], file(inputFile))
    }

    NORMALISE_VCF(channel.fromList(vcfInputs))

    def byTool = NORMALISE_VCF.out.vcf.flatMap { chMeta, vcf, tbi, excluded ->
        def chConf = challenges[chMeta.challenge]
        chConf.tools.findAll { _name, toolConf -> toolConf.enabled }.collect { toolName, toolConf ->
            tuple(toolName, toolMeta(chMeta.challenge, chConf, toolName, toolConf), vcf, tbi, excluded)
        }
    }

    KING_PLINK(forTool(byTool, 'king_plink'))
    NGSRELATE(forTool(byTool, 'ngsrelate'))
    ADMIXTURE(forTool(byTool, 'admixture'))
    SHAPEIT5(forTool(byTool, 'shapeit5'))
    READ_KIN(forTool(byTool, 'read_kin'))
    ICMLKIN(forTool(byTool, 'icmlkin'))
    GRUPS_RS(forTool(byTool, 'grups_rs'))

    emit:
    king_plink = KING_PLINK.out.result
    ngsrelate  = NGSRELATE.out.result
    admixture  = ADMIXTURE.out.result
    shapeit5   = SHAPEIT5.out.result
}
