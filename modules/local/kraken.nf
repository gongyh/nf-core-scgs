process KRAKEN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kraken2=2.1.2 bioconda::krona=2.7.1 bioconda::krakentools=1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-3bbb1b9ff2130265cf8d9498a097b04978fb988f:6688dcb6662e35001e709b425821fff321f15540-0' :
        'scgs/mulled-v2-3bbb1b9ff2130265cf8d9498a097b04978fb988f:6688dcb6662e35001e709b425821fff321f15540-0' }"

    input:
    tuple val(meta), path(reads)
    path db
    path taxonomy, stageAs: 'taxonomy.tab'

    output:
    tuple val(meta), path("*.krk"),    emit: report
    tuple val(meta), path("*.html"),   emit: html
    path "versions.yml",               emit: versions

    script:
    def mode = meta.single_end ? "" : "--paired"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    TAXONOMY=\$(find -L . -name '*.tab' -exec dirname {} \\;)
    kraken2 --db $db --threads ${task.cpus} --report ${prefix}.krk --output ${prefix}.k2 --gzip-compressed ${mode} $reads
    kreport2krona.py -r ${prefix}.krk -o ${prefix}.krn
    ktImportText -o ${prefix}_taxonomy.html ${prefix}.krn

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken --version 2>&1) | sed 's/^.*kraken //; s/Using.*\$//')
    END_VERSIONS
    """
}
