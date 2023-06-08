process KRAKEN {
    tag "$meta.id"
    publishDir "${params.outdir}/kraken", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.report")     , emit: report
    tuple val(meta), path("*.krona.html") , emit: html

    when:
    params.kraken_db

    script:
    def mode = meta.single_end ? "" : "--paired"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kraken -db $db --threads ${task.cpus} --fastq-input --gzip-compressed ${mode} --check-names --output ${prefix}.krk $reads
    kraken-report -db $db ${prefix}.krk > ${prefix}.report
    cut -f2,3 ${prefix}.krk > ${prefix}.f23
    ktImportTaxonomy -o ${prefix}.krona.html ${prefix}.f23
    """
}
