process METARON {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs)
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("$prefix"), emit: out_put
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    metaron.py -n ${prefix} -p op -i ${tab} -j ${contigs} -t 2 -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
