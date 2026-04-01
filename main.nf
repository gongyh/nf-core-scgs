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

include { SCGS } from './workflows/scgs'
include { MINIMETA } from './workflows/minimeta'

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW BASED ON PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    if (params.minimeta) {
        NFCORE_MINIMETA ()
    } else {
        NFCORE_SCGS ()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
