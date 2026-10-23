/*
 * Sweep Detection category workflow. Same fan-out pattern as
 * workflows/genetic_history.nf - see that file's header comment. Every
 * sweep challenge is VCF-based and passes through NORMALISE_VCF first.
 */

include { NORMALISE_VCF } from '../modules/common/normalise_vcf/main.nf'
include { SWEEPFINDER2  } from '../modules/sweep_detection/sweepfinder2/main.nf'
include { FLEXSWEEP     } from '../modules/sweep_detection/flexsweep/main.nf'
include { SWEED         } from '../modules/sweep_detection/sweed/main.nf'
include { OMEGAPLUS     } from '../modules/sweep_detection/omegaplus/main.nf'
include { RAISD         } from '../modules/sweep_detection/raisd/main.nf'
include { TIMESWEEPER   } from '../modules/sweep_detection/timesweeper/main.nf'

def toolMeta(String chName, Map chConf, String _toolName, Map toolConf) {
    return [
        category      : 'sweep_detection',
        challenge     : chName,
        accessible_bed: chConf.accessible_bed,
        conda         : toolConf.conda,
        args          : toolConf.args ?: '',
        options       : toolConf.options ?: [:],
        reason        : toolConf.reason ?: '',
    ]
}

// Entries of the (tool, meta, vcf, tbi, excluded_bed) fan-out channel for one tool.
def forTool(byTool, String name) {
    return byTool.filter { entry -> entry[0] == name }
}

// SweepFinder2/FlexSweep additionally take the challenge's accessible-region mask.
def withAccessible(byTool, String name) {
    return forTool(byTool, name).map { _tool, meta, vcf, _tbi, _excluded ->
        tuple(meta, vcf, file(meta.accessible_bed, checkIfExists: true))
    }
}

def vcfOnly(byTool, String name) {
    return forTool(byTool, name).map { _tool, meta, vcf, _tbi, _excluded -> tuple(meta, vcf) }
}

workflow SWEEP_DETECTION {
    main:
    def challenges = params.sweep_detection.challenges

    def vcfInputs = []
    challenges.each { chName, chConf ->
        if (!chConf.enabled) return
        if (!chConf.tools.any { _name, toolConf -> toolConf.enabled }) return
        def inputFile = GhistUtils.testingFile(params.data_dir, chConf.testing_vcf_name)
        if (!inputFile.exists()) {
            log.warn "sweep_detection/${chName}: testing file not found (${inputFile}), skipping challenge"
            return
        }
        vcfInputs << tuple([category: 'sweep_detection', challenge: chName], file(inputFile))
    }

    NORMALISE_VCF(channel.fromList(vcfInputs))

    def byTool = NORMALISE_VCF.out.vcf.flatMap { chMeta, vcf, tbi, excluded ->
        def chConf = challenges[chMeta.challenge]
        chConf.tools.findAll { _name, toolConf -> toolConf.enabled }.collect { toolName, toolConf ->
            tuple(toolName, toolMeta(chMeta.challenge, chConf, toolName, toolConf), vcf, tbi, excluded)
        }
    }

    SWEEPFINDER2(withAccessible(byTool, 'sweepfinder2'))
    FLEXSWEEP(withAccessible(byTool, 'flexsweep'))
    SWEED(vcfOnly(byTool, 'sweed'))
    OMEGAPLUS(vcfOnly(byTool, 'omegaplus'))
    RAISD(vcfOnly(byTool, 'raisd'))
    TIMESWEEPER(vcfOnly(byTool, 'timesweeper'))

    emit:
    sweepfinder2 = SWEEPFINDER2.out.result
    flexsweep    = FLEXSWEEP.out.result
}
