process KRAKEN {
    tag "$prefix"
    publishDir "${params.outdir}/kraken", mode: 'copy'

    input:
    path(reads)
    path db
    val(single_end)

    output:
    path("${prefix}.report"),                            emit: report
    path("${prefix}.krona.html")

    when:
    params.kraken_db

    script:
    prefix = reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    def mode = single_end ? "" : "--paired"
    """
    kraken -db $db --threads ${task.cpus} --fastq-input --gzip-compressed ${mode} --check-names --output ${prefix}.krk $reads
    kraken-report -db $db ${prefix}.krk > ${prefix}.report
    cut -f2,3 ${prefix}.krk > ${prefix}.f23
    ktImportTaxonomy -o ${prefix}.krona.html ${prefix}.f23
    """
}
