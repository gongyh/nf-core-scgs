process MULTIQC_REF {
    label "multiqc"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path multiqc_config
    path('fastqc/*')
    path('trimgalore/*')
    path('software_versions/*')
    path('samtools/*')
    path('preseq/*')
    path('*')
    path('quast/*')
    path('prokka/*')
    path('kraken/*')
    path workflow_summary

    output:
    path("*multiqc_report.html"), emit: multiqc_report1
    path("*_data")

    when:
    denovo == false

    script:
    """
    multiqc -f --config $multiqc_config .
    """
}
