process KRAKEN {
    tag "$meta.id"

    conda "kraken=1.1.1,krona=2.7.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d45752891fea2584428a164c55ff535957eb7fa2:17bc7e8d082e77491b01a53af02d08779b923f10-0' :
        'scgs/mulled-v2-d45752891fea2584428a164c55ff535957eb7fa2:17bc7e8d082e77491b01a53af02d08779b923f10-0' }"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.report")     , emit: report
    tuple val(meta), path("*.krona.html") , emit: html
    path "versions.yml",                    emit: versions

    script:
    def mode = meta.single_end ? "" : "--paired"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kraken -db $db --threads ${task.cpus} --fastq-input --gzip-compressed ${mode} --check-names --output ${prefix}.krk $reads
    kraken-report -db $db ${prefix}.krk > ${prefix}.report
    cut -f2,3 ${prefix}.krk > ${prefix}.f23
    ktImportTaxonomy -o ${prefix}.krona.html ${prefix}.f23

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken: \$(echo \$(kraken --version 2>&1) | sed 's/^.*kraken //; s/Using.*\$//')
    END_VERSIONS
    """
}
