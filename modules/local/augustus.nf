process AUGUSTUS {
    tag "$prefix"
    publishDir "${params.outdir}/augustus", mode: 'copy'

    input:
    path contigs

    output:
    path("${prefix}.aa"),                      emit: faa_eukcc
    path("${prefix}*")

    when:
    euk

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    """
    # clean id
    cat $contigs | sed 's/_length.*\$//g' > ${prefix}_clean.fasta
    # mask genome
    tantan ${prefix}_clean.fasta > ${prefix}_mask.fasta
    # gene prediction
    augustus --species=${params.augustus_species} --gff3=on --uniqueGeneId=true --protein=on --codingseq=on ${prefix}_mask.fasta > ${prefix}.gff
    # generate proteins
    getAnnoFasta.pl ${prefix}.gff
    """
}
