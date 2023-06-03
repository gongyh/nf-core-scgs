process SPADES {
    tag "${prefix}"
    publishDir "${params.outdir}/spades", mode: 'copy'

    input:
    path clean_reads
    val(single_end)

    output:
    path("${prefix}.ctg200.fasta"),                    emit: ctg200
    path("${prefix}.ctgs.fasta"),                      emit: ctg

    when:
    params.ass

    script:
    prefix = clean_reads[0].toString() - ~/(_trimmed)?(_norm)?(_combined)?(\.R1)?(_1)?(_R1)?(\.1_val_1)?(_1_val_1)?(_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?(\.bz2)?$/
    R1 = clean_reads[0].toString()
    def mode = params.bulk ? "bulk" : "mda"
    if (single_end) {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
    spades.py -s $R1 --careful --cov-cutoff auto -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    else
    spades.py --sc -s $R1 --careful -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    fi
    ln -s ${prefix}.spades_out/contigs.fasta ${prefix}.contigs.fasta
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/_length.*\$//g' > ${prefix}.ctgs.fasta
    """
    } else {
    R2 = clean_reads[1].toString()
    """
    if [ \"${mode}\" == \"bulk\" ]; then
    spades.py -1 $R1 -2 $R2 --careful --cov-cutoff auto -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    else
    spades.py --sc -1 $R1 -2 $R2 --careful -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    fi
    ln -s ${prefix}.spades_out/contigs.fasta ${prefix}.contigs.fasta
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/_length.*\$//g' > ${prefix}.ctgs.fasta
    """
    }
}
