process TSNE {
    tag "$meta.id"
    publishDir "${params.outdir}/tsne", mode: 'copy'

    input:
    tuple val(meta), path(contigs)

    output:
    path("${prefix}_tsne.tsv")

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
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
