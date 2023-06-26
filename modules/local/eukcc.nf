process EUKCC {
    tag "$meta.id"
    label 'process_medium'

    conda "eukcc=2.1.0--pypyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukcc:2.1.0--pypyhdfd78af_0' :
        'biocontainers/eukcc:2.1.0--pypyhdfd78af_0' }"

    input:
    tuple val(meta), path(contig)
    path db

    output:
    path("${prefix}")

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    eukcc single --out $prefix --db $db --threads 8 $contig || echo "Ignore minor errors of eukcc!"
    """
}
