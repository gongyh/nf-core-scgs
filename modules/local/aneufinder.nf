process ANEUFINDER {
    publishDir "${pp_outdir}", mode: 'copy'

    input:
    path("bams/*")
    path("bams/*")

    output:
    path('CNV_output'),  emit: cnv

    when:
    !params.bulk && params.cnv && !single_end && !params.nanopore

    script:
    pp_outdir = "${params.outdir}/aneufinder"
    """
    aneuf.R ./bams CNV_output ${task.cpus}
    """
}
