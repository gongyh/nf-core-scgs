process KTUPDATETAXONOMY {
    label 'process_single'

    conda "bioconda::kraken2=2.1.2 bioconda::krona=2.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b85de0f0888e1a8481d8c5d0c3b52736036932cc:96c1f81ca967332ad179c5bc0a350133f6bdf2a8-0' :
        'scgs/mulled-v2-b85de0f0888e1a8481d8c5d0c3b52736036932cc:96c1f81ca967332ad179c5bc0a350133f6bdf2a8-0' }"

    output:
    path 'taxonomy/taxonomy.tab', emit: taxonomy
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ktUpdateTaxonomy.sh \\
        $args \\
        taxonomy/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: '2.7.1'
    END_VERSIONS
    """
}
