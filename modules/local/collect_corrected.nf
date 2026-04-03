process COLLECT_CORRECTED {
    tag "collect"
    label 'process_medium'
    publishDir "${params.outdir}/joint", mode: 'copy'

    input:
    path p1_files   // 所有子样本的 P1_corrected 文件列表
    path p2_files   // 所有子样本的 P2_corrected 文件列表
    path s_files    // 所有子样本的 S_corrected 文件列表

    output:
    path "all_corrected_R1.fastq.gz", emit: r1
    path "all_corrected_R2.fastq.gz", emit: r2
    path "all_corrected_S.fastq.gz",  emit: s
    path "versions.yml", emit: versions

    script:
    """
    cat ${p1_files.join(' ')} > all_corrected_R1.fastq.gz
    cat ${p2_files.join(' ')} > all_corrected_R2.fastq.gz
    cat ${s_files.join(' ')}  > all_corrected_S.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        script: "collect_corrected.nf"
    END_VERSIONS
    """
}
