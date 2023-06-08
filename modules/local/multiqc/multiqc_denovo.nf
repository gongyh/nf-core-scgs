process MULTIQC_DENOVO {
    label "multiqc"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path multiqc_config
    path('fastqc/*')
    path('trimgalore/*')
    path('software_versions/*')
    path('quast/*')
    path('prokka/*')
    path('kraken/*')
    path workflow_summary

    output:
    path("*multiqc_report.html"), emit: multiqc_report2
    path("*_data")

    when:
    denovo == true

    script:
    """
    multiqc -f --config $multiqc_config .
    """
}
