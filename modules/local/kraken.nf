process KRAKEN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kraken2=2.1.2 bioconda::krona=2.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b85de0f0888e1a8481d8c5d0c3b52736036932cc:96c1f81ca967332ad179c5bc0a350133f6bdf2a8-0' :
        'scgs/mulled-v2-b85de0f0888e1a8481d8c5d0c3b52736036932cc:96c1f81ca967332ad179c5bc0a350133f6bdf2a8-0' }"

    input:
    tuple val(meta), path(reads)
    path db
    path taxonomy, stageAs: 'taxonomy.tab'

    output:
    tuple val(meta), path("*.report")     , emit: report
    tuple val(meta), path("*.krona.html") , emit: html
    path "versions.yml",                    emit: versions

    script:
    def mode = meta.single_end ? "" : "--paired"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    TAXONOMY=\$(find -L . -name '*.tab' -exec dirname {} \\;)
    kraken2 --db $db --threads ${task.cpus} --report ${prefix}.report --output ${prefix}.krk --gzip-compressed ${mode} --use-names $reads
    cut -f2,3 ${prefix}.krk > ${prefix}.f23
    ktImportTaxonomy -o ${prefix}.krona.html -tax \$TAXONOMY/ ${prefix}.f23

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken --version 2>&1) | sed 's/^.*kraken //; s/Using.*\$//')
    END_VERSIONS
    """
}
