process ANEUFINDER {
    label 'process_medium'

    conda "bioconda::bioconductor-aneufinder=1.26.0=r42hf17093f_1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aneufinder:1.26.0--r42hf17093f_1' :
        'biocontainers/bioconductor-aneufinder:1.26.0--r42hf17093f_1' }"

    input:
    path("bams/*")
    path("bams/*")

    output:
    path('CNV_output') , emit: cnv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    aneuf.R ./bams CNV_output ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aneufinder: \$(Rscript -e 'v=format(packageVersion("AneuFinder"));cat(v)')
    END_VERSIONS
    """
}
