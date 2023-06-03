process REMAP {
    tag "${prefix}"
    publishDir "${params.outdir}/remap", mode: 'copy'

    input:
    path reads
    path index
    val(single_end)

    output:
    path("${prefix}_ass.bam"),                         emit: bam

    when:
    params.remap

    script:
    prefix = reads[0].toString() - ~/(_trimmed)?(_norm)?(_combined)?(\.R1)?(_1)?(_R1)?(\.1_val_1)?(_1_val_1)?(_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    R1 = reads[0].toString()
    def filtering = params.allow_multi_align ? '' : "| samtools view -b -q 40 -F 4 -F 256 -"
    if (single_end) {
    """
    bowtie2 -x ${prefix}Bowtie2Index/${prefix} -p ${task.cpus} -U $R1 | samtools view -bT ${prefix}Bowtie2Index - $filtering > ${prefix}_ass.bam
    """
    } else {
    R2 = reads[1].toString()
    """
    bowtie2 --no-mixed --no-discordant -X 1000 -x ${prefix}Bowtie2Index/${prefix} -p ${task.cpus} -1 $R1 -2 $R2 | samtools view -bT ${prefix}Bowtie2Index - $filtering > ${prefix}_ass.bam
    """
    }
}
