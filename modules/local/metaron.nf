process METARON {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs)
    path(gene_models)

    output:
    tuple val(meta), path("$prefix"), emit: out_operon
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    metaron.py -n ${prefix} -p op -i ${prefix}.gff -j ${contigs} -t 2 -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaRon: 1.0
    END_VERSIONS
    """
}
