process REMAP {
    tag "$meta.id"
    label 'process_medium'

    /**
    conda "bioconda::bowtie2=2.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0' :
        'biocontainers/bowtie2:2.4.4--py39hbb4e92a_0' }"
    */

    input:
    tuple val(meta), path(reads)
    path(index)
    val(allow_multi_align)

    output:
    tuple val(meta), path("${prefix}_ass.bam"),  emit: bam

    when:
    params.remap

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def filtering = allow_multi_align ? '' : "| samtools view -b -q 40 -F 4 -F 256 -"
    if (meta.single_end) {
    """
    bowtie2 -x ${prefix}Bowtie2Index/${prefix} -p ${task.cpus} -U ${reads[0]} | samtools view -bT ${prefix}Bowtie2Index - $filtering > ${prefix}_ass.bam
    """
    } else {
    """
    bowtie2 --no-mixed --no-discordant -X 1000 -x ${prefix}Bowtie2Index/${prefix} -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} | samtools view -bT ${prefix}Bowtie2Index - $filtering > ${prefix}_ass.bam
    """
    }
}
