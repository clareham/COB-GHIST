/*
 * Sweep Detection category workflow. Same fan-out pattern as
 * workflows/genetic_history.nf - see that file's header comment.
 */

include { SWEEPFINDER2 } from '../modules/sweep_detection/sweepfinder2/main.nf'
include { FLEXSWEEP    } from '../modules/sweep_detection/flexsweep/main.nf'
include { SWEED        } from '../modules/sweep_detection/sweed/main.nf'
include { OMEGAPLUS    } from '../modules/sweep_detection/omegaplus/main.nf'
include { RAISD        } from '../modules/sweep_detection/raisd/main.nf'
include { TIMESWEEPER  } from '../modules/sweep_detection/timesweeper/main.nf'

workflow SWEEP_DETECTION {
    main:
    def byTool = [:].withDefault { [] }

    params.sweep_detection.challenges.each { chName, chConf ->
        if (!chConf.enabled) return
        def inputFile = GhistUtils.testingFile(params.data_dir, chConf.testing_vcf_name)
        if (!inputFile.exists()) {
            log.warn "sweep_detection/${chName}: testing file not found (${inputFile}), skipping challenge"
            return
        }
        chConf.tools.each { toolName, toolConf ->
            if (!toolConf.enabled) return
            def meta = [
                category      : 'sweep_detection',
                challenge     : chName,
                accessible_bed: chConf.accessible_bed,
                conda         : toolConf.conda,
                args          : toolConf.args ?: '',
                options       : toolConf.options ?: [:],
                reason        : toolConf.reason ?: '',
            ]
            if (toolName in ['sweepfinder2', 'flexsweep']) {
                def accessible = file(chConf.accessible_bed, checkIfExists: true)
                byTool[toolName] << tuple(meta, file(inputFile), accessible)
            } else {
                byTool[toolName] << tuple(meta, file(inputFile))
            }
        }
    }

    SWEEPFINDER2(channel.fromList(byTool['sweepfinder2']))
    FLEXSWEEP(channel.fromList(byTool['flexsweep']))
    SWEED(channel.fromList(byTool['sweed']))
    OMEGAPLUS(channel.fromList(byTool['omegaplus']))
    RAISD(channel.fromList(byTool['raisd']))
    TIMESWEEPER(channel.fromList(byTool['timesweeper']))

    emit:
    sweepfinder2 = SWEEPFINDER2.out.result
    flexsweep    = FLEXSWEEP.out.result
}
