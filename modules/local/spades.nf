process SPADES {
    tag "${meta.id}"
    publishDir "${params.outdir}/spades", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.ctg200.fasta"),   emit: ctg200
    tuple val(meta), path("${prefix}.ctgs.fasta"),     emit: ctg

    when:
    params.ass

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def mode = params.bulk ? "bulk" : "mda"
    if (meta.single_end) {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
    spades.py -s ${reads[0]} --careful --cov-cutoff auto -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    else
    spades.py --sc -s ${reads[0]} --careful -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    fi
    ln -s ${prefix}.spades_out/contigs.fasta ${prefix}.contigs.fasta
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/_length.*\$//g' > ${prefix}.ctgs.fasta
    """
    } else {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
    spades.py -1 ${reads[0]} -2 ${reads[1]} --careful --cov-cutoff auto -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    else
    spades.py --sc -1 ${reads[0]} -2 ${reads[1]} --careful -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    fi
    ln -s ${prefix}.spades_out/contigs.fasta ${prefix}.contigs.fasta
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/_length.*\$//g' > ${prefix}.ctgs.fasta
    """
    }
}