process REMAP {
    tag "${prefix}"
    publishDir "${params.outdir}/remap", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(index)

    output:
    tuple val(meta), path("${prefix}_ass.bam"),  emit: bam

    when:
    params.remap

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def filtering = params.allow_multi_align ? '' : "| samtools view -b -q 40 -F 4 -F 256 -"
    if (single_end) {
    """
    bowtie2 -x ${prefix}Bowtie2Index/${prefix} -p ${task.cpus} -U ${reads[0]} | samtools view -bT ${prefix}Bowtie2Index - $filtering > ${prefix}_ass.bam
    """
    } else {
    """
    bowtie2 --no-mixed --no-discordant -X 1000 -x ${prefix}Bowtie2Index/${prefix} -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} | samtools view -bT ${prefix}Bowtie2Index - $filtering > ${prefix}_ass.bam
    """
    }
}