process MULTIQC_DENOVO {
    label "multiqc"

    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"

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
    path("*multiqc_report.html"), emit: report
    path("*_data")

    script:
    """
    multiqc -f --config $multiqc_config .
    """
}
