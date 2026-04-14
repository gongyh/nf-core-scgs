process REMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bowtie2=2.4.4 bioconda::samtools=1.17"
    container "scgs/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:7a723346025146e191fdbf519b8131c258b9eeab-0"

    input:
    tuple val(meta), path(reads), path(index)
    val(allow_multi_align)

    output:
    tuple val(meta), path("${prefix}_ass.sort.bam")    , emit: bam
    tuple val(meta), path("${prefix}_ass.sort.bam.bai"), emit: bai
    tuple val(meta), path("${prefix}_ass.sort.bam"), path("${prefix}_ass.sort.bam.bai"), emit: bam_bai
    path "versions.yml"                                , emit: versions

    when:true
    println "REMAP process started for ${meta.id}"
    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    index_prefix = index.find { it.toString().endsWith(".1.bt2") }?.toString()?.replaceAll("\\.1.bt2\$", "")
    if (!index_prefix) error "Index files not found in ${index}"
    def filter_args = allow_multi_align ? '' : '| samtools view -b -q 40 -F 4 -F 256 -'
    if (meta.single_end) {
    """
    bowtie2 -x ${index_prefix} -p ${task.cpus} -U ${reads} | samtools view -bS ${filter_args} > ${prefix}_ass.bam
    samtools sort -o ${prefix}_ass.sort.bam ${prefix}_ass.bam
    samtools index ${prefix}_ass.sort.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    } else {
    """
    bowtie2 --no-mixed --no-discordant -X 1000 -x ${index_prefix} -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} | samtools view -bS - ${filter_args} > ${prefix}_ass.bam
    samtools sort -o ${prefix}_ass.sort.bam ${prefix}_ass.bam
    samtools index ${prefix}_ass.sort.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    }
}
