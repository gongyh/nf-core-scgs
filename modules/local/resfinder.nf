process RESFINDER {
    tag "$prefix"
    publishDir "${params.outdir}/ARG", mode: 'copy'

    input:
    path contigs
    path db

    output:
    path("${prefix}/*")

    when:
    !euk && params.acquired && resfinder_db

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    """
    mkdir -p $prefix
    python /opt/resfinder/resfinder.py -i $contigs -o $prefix -p $db -mp blastn -x
    rm -rf $prefix/tmp
    """
}
