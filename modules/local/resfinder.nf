process RESFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::resfinder==4.1.11--hdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/resfinder:4.1.11--hdfd78af_0' :
        'biocontainers/resfinder:4.1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(contigs)
    path db

    output:
    path("${prefix}/*")

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $prefix
    python /opt/resfinder/resfinder.py -i $contigs -o $prefix -p $db -mp blastn -x
    rm -rf $prefix/tmp
    """
}
