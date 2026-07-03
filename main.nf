#!/usr/bin/env nextflow

/*
========================================================================================
                        gongyh/nf-core-scgs
========================================================================================
    gongyh/nf-core-scgs Analysis Pipeline.
    #### Homepage / Documentation
    https://github.com/gongyh/nf-core-scgs
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SCGS              } from './workflows/scgs'
include { MINIMETA          } from './workflows/minimeta'
include { PREPARE_DATABASES } from './workflows/prepare_databases'

include { completionEmail   } from './subworkflows/nf-core/utils_nfcore_pipeline/main'
include { completionSummary } from './subworkflows/nf-core/utils_nfcore_pipeline/main'

include { helpMessage as helpMessageSCGS           } from './workflows/scgs'
include { helpMessage as helpMessageMinimeta       } from './workflows/minimeta'
include { helpMessage as helpMessagePrepareDB      } from './workflows/prepare_databases'

//
// WORKFLOW: Run SCGS analysis pipeline
//

workflow NFCORE_SCGS {
    SCGS ()
}

//
// WORKFLOW: Run MINIMETA analysis pipeline
//

workflow NFCORE_MINIMETA {
    MINIMETA ()
}

//
// WORKFLOW: Run database preparation pipeline
//

workflow NFCORE_PREPARE_DATABASES {
    PREPARE_DATABASES ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW BASED ON PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    if (params.prepare_databases) {
        // Show help message
        if (params.help){
            helpMessagePrepareDB()
            exit 0
        }
        NFCORE_PREPARE_DATABASES ()
    } else if (params.minimeta) {
        // Show help message
        if (params.help){
            helpMessageMinimeta()
            exit 0
        }
        NFCORE_MINIMETA ()
    } else {
        // Show help message
        if (params.help){
            helpMessageSCGS()
            exit 0
        }
        NFCORE_SCGS ()
    }
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {
    if (params.email){
        completionEmail(summary_params,
            params.email,
            null,
            false,
            params.outdir,
            log,
            multiqc_report.getVal()
        )
    }
    completionSummary()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
