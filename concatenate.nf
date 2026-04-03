process CONCATENATE {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/import", mode: 'copy'

    input:
    tuple val(meta), path(reads)   // 双端时 reads 为 [r1, r2]

    output:
    tuple val(meta), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz"), emit: merged_reads
    path "versions.yml", emit: versions

    script:
    if (meta.single_end) {
        """
        ln -sf ${reads} ${meta.id}_R1.fastq.gz
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            script: "concatenate.nf (passthrough)"
        END_VERSIONS
        """
    } else {
        """
        ln -sf ${reads[0]} ${meta.id}_R1.fastq.gz
        ln -sf ${reads[1]} ${meta.id}_R2.fastq.gz
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            script: "concatenate.nf (passthrough)"
        END_VERSIONS
        """
    }
}
