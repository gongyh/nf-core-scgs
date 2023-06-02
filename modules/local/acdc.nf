process ACDC {
    tag "${prefix}"
    publishDir "${params.outdir}/acdc", mode: 'copy'

    input:
    path contigs
    path db
    path tax

    output:
    path("${prefix}")

    when:
    false

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    """
    cat $tax | grep -v '^#' | cut -f1,18 > genus.txt
    /usr/local/bin/acdc -i $contigs -m 1000 -b 100 -o $prefix -K $db -x genus.txt -T ${task.cpus}
    """
}
