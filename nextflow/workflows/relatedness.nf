/*
 * Relatedness Estimation category workflow. Same fan-out pattern as
 * workflows/genetic_history.nf - see that file's header comment.
 *
 * relatedness1 contributes no entries here: its challenge config has an
 * empty `tools {}` map (data/tool format mismatch - see
 * conf/relatedness/relatedness1.config).
 */

include { KING_PLINK } from '../modules/relatedness/king_plink/main.nf'
include { NGSRELATE  } from '../modules/relatedness/ngsrelate/main.nf'
include { ADMIXTURE  } from '../modules/relatedness/admixture/main.nf'
include { SHAPEIT5   } from '../modules/relatedness/shapeit5/main.nf'
include { READ_KIN   } from '../modules/relatedness/read_kin/main.nf'
include { ICMLKIN    } from '../modules/relatedness/icmlkin/main.nf'
include { GRUPS_RS   } from '../modules/relatedness/grups_rs/main.nf'

workflow RELATEDNESS {
    main:
    def byTool = [:].withDefault { [] }

    params.relatedness.challenges.each { chName, chConf ->
        if (!chConf.enabled) return
        def fileName = chConf.input_type == 'vcf' ? chConf.testing_vcf_name : chConf.testing_csv_name
        def inputFile = GhistUtils.testingFile(params.data_dir, fileName)
        if (!inputFile.exists()) {
            log.warn "relatedness/${chName}: testing file not found (${inputFile}), skipping challenge"
            return
        }
        chConf.tools.each { toolName, toolConf ->
            if (!toolConf.enabled) return
            if (toolName in ['king_plink', 'admixture'] &&
                (!(chConf.autosomes instanceof Integer) || chConf.autosomes < 1)) {
                error "relatedness/${chName}: ${toolName} requires a positive diploid autosome count"
            }
            def meta = [
                category  : 'relatedness',
                challenge : chName,
                input_type: chConf.input_type,
                autosomes : chConf.autosomes,
                conda     : toolConf.conda,
                args      : toolConf.args ?: '',
                k         : toolConf.k ?: 3,
                reason    : toolConf.reason ?: '',
            ]
            byTool[toolName] << tuple(meta, file(inputFile))
        }
    }

    KING_PLINK(channel.fromList(byTool['king_plink']))
    NGSRELATE(channel.fromList(byTool['ngsrelate']))
    ADMIXTURE(channel.fromList(byTool['admixture']))
    SHAPEIT5(channel.fromList(byTool['shapeit5']))
    READ_KIN(channel.fromList(byTool['read_kin']))
    ICMLKIN(channel.fromList(byTool['icmlkin']))
    GRUPS_RS(channel.fromList(byTool['grups_rs']))

    emit:
    king_plink = KING_PLINK.out.result
    ngsrelate  = NGSRELATE.out.result
    admixture  = ADMIXTURE.out.result
    shapeit5   = SHAPEIT5.out.result
}
