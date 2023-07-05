process ACDC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::acdc=1.02=h4ac6f70_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/acdc:1.02--h4ac6f70_0' :
        'biocontainers/acdc:1.02--h4ac6f70_0' }"

    input:
    tuple val(meta), path(contigs)
    tuple val(meta), path(tax)
    path db

    output:
    tuple val(meta), path("${prefix}") , emit: out_put
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $tax | grep -v '^#' | cut -f1,18 > genus.txt
    /usr/local/bin/acdc -i $contigs -m 1000 -b 100 -o $prefix -K $db -x genus.txt -T ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        acdc: '1.02'
    END_VERSIONS
    """
}
