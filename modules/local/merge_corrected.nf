process MERGE_CORRECTED {
    tag "merge"
    label 'process_low'
    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    path p1_files
    path p2_files

    output:
    path "all_R1.fastq.gz", emit: r1
    path "all_R2.fastq.gz", emit: r2
    path "versions.yml", emit: versions

    script:
    """
    cat ${p1_files.join(' ')} > all_R1.fastq.gz
    cat ${p2_files.join(' ')} > all_R2.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        script: "merge_corrected.nf"
    END_VERSIONS
    """
}
