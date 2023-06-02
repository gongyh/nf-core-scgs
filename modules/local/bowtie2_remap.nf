process BOWTIE2_REMAP {
    tag "${prefix}"
    publishDir "${params.outdir}/remap_bowtie2_index", mode: 'copy'

    input:
    path contigs

    output:
    path("${prefix}Bowtie2Index"),                 emit: index

    when:
    params.remap

    script:
    prefix = contigs.toString() - ~/(\.ctg200\.fasta?)(\.ctgs\.fasta)?(\.ctgs)?(\.ctg200)?(\.fasta)?(\.fa)?$/
    """
    mkdir -p ${prefix}Bowtie2Index; cd ${prefix}Bowtie2Index
    ln -s ../${contigs} ${prefix}.fa
    bowtie2-build ${prefix}.fa ${prefix}
    """
}
