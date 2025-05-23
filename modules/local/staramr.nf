process STARAMR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda:staramr=0.8.0=pyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staramr:0.8.0--pyhdfd78af_0' :
        'biocontainers/staramr:0.8.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs)
    val(acquired)
    val(point)
    val(species)

    output:
    tuple val(meta), path("${prefix}/*"), emit: out_put
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def species = species
    if (acquired && !point) {
    """
    staramr search -o $prefix $contigs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staramr: \$(echo \$(staramr -V 2>&1) | sed 's/^.*staramr //; s/Using.*\$//')
    END_VERSIONS
    """
    } else {
    """
    staramr search --pointfinder-organism $species -o $prefix $contigs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staramr: \$(echo \$(staramr -V 2>&1) | sed 's/^.*staramr //; s/Using.*\$//')
    END_VERSIONS
    """
    }
}
