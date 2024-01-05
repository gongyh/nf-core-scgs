process GRAPHBIN {
    label 'process_medium'

    conda "bioconda::graphbin=1.7.1--pyh7cba7a3_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphbin:1.7.1--pyh7cba7a3_0' :
        'biocontainers/graphbin:1.7.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(contig)
    tuple val(meta), path(path)
    tuple val(meta), path(gfa)
    path("*")

    output:
    path("binning/*")      , emit: out_put
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p binning/${meta.id}
    graphbin --assembler spades --graph $gfa --contigs $contig --paths $path --binned ${prefix}.bin.csv --output binning/${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphbin: \$(echo \$(graphbin -v 2>&1) | sed 's/^.*graphbin, version //; s/Using.*\$//')
    END_VERSIONS
    """
}
