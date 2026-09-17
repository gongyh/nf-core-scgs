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
    main:
    SCGS ()

    emit:
    summary_params = SCGS.out.summary_params
    multiqc_report = SCGS.out.multiqc_report
    published = SCGS.out.published
}

//
// WORKFLOW: Run MINIMETA analysis pipeline
//

workflow NFCORE_MINIMETA {
    main:
    MINIMETA ()

    emit:
    summary_params = MINIMETA.out.summary_params
    multiqc_report = MINIMETA.out.multiqc_report
    published = MINIMETA.out.published
}

//
// WORKFLOW: Run database preparation pipeline
//

workflow NFCORE_PREPARE_DATABASES {
    main:
    PREPARE_DATABASES ()

    emit:
    published = PREPARE_DATABASES.out.published
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW BASED ON PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    main:
    ch_published = channel.empty()

    if (params.prepare_databases) {
        // Show help message
        if (params.help){
            helpMessagePrepareDB()
            exit 0
        }
        NFCORE_PREPARE_DATABASES ()
        ch_published = NFCORE_PREPARE_DATABASES.out.published
        workflow.onComplete = {
            completionSummary()
        }
    } else if (params.minimeta) {
        // Show help message
        if (params.help){
            helpMessageMinimeta()
            exit 0
        }
        NFCORE_MINIMETA ()
        ch_published = NFCORE_MINIMETA.out.published
        workflow.onComplete = {
            if (params.email) {
                completionEmail(
                    NFCORE_MINIMETA.out.summary_params.getVal(),
                    params.email,
                    null,
                    false,
                    params.outdir,
                    log,
                    NFCORE_MINIMETA.out.multiqc_report.getVal()
                )
            }
            completionSummary()
        }
    } else {
        // Show help message
        if (params.help){
            helpMessageSCGS()
            exit 0
        }
        NFCORE_SCGS ()
        ch_published = NFCORE_SCGS.out.published
        workflow.onComplete = {
            if (params.email) {
                completionEmail(
                    NFCORE_SCGS.out.summary_params.getVal(),
                    params.email,
                    null,
                    false,
                    params.outdir,
                    log,
                    NFCORE_SCGS.out.multiqc_report.getVal()
                )
            }
            completionSummary()
        }
    }

    publish:
    pipeline_results = ch_published
}

output {
    pipeline_results {
        path { result -> "${result.destination}/" }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
