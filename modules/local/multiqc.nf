process MULTIQC {
    label 'process_medium'

    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    path multiqc_config
    path('fastqc/*')
    path('trimgalore/*')
    path('software_versions/*')
    path('fastqc2/*')
    path('samtools/*')
    path('preseq/*')
    path('*')
    path('quast/*')
    path('checkm/*')
    path('prokka/*')
    path('kraken/*')
    path workflow_summary

    output:
    path("*multiqc_report.html"), emit: report
    path("*_data")

    script:
    """
    multiqc -f --config $multiqc_config .
    """
}