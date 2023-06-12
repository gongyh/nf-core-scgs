process TSNE {
    tag "${prefix}"
    publishDir "${params.outdir}/tsne", mode: 'copy'

    input:
    path contigs

    output:
    path("${prefix}_tsne.tsv")

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    """
    faFilterByLen.pl ${contigs} 1000 > ${prefix}.ctg1k.fasta
    if [ -s ${prefix}.ctg1k.fasta ]
    then
    kpal count -k 4 -r ${prefix}.ctg1k.fasta ${prefix}.4mer
    kmer_tsne.py ${prefix}.4mer ${prefix}_tsne.tsv ${task.cpus}
    else
    touch ${prefix}_tsne.tsv
    fi
    """
}
