nextflow.enable.types = true

process STARAMR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda:staramr=0.8.0=pyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staramr:0.8.0--pyhdfd78af_0' :
        'biocontainers/staramr:0.8.0--pyhdfd78af_0' }"

    input:
    tuple(meta: Map, contigs: Path, acquired: Boolean, point: Boolean, species: String)

    output:
    record(meta: meta, out_put: file("${task.ext.prefix ?: meta.id}", type: "dir"), versions: file('versions.yml'))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
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
