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
    path "manifest.txt", emit: manifest
    path "versions.yml", emit: versions

    script:
    """
    p1_sorted=\$(printf '%s\\n' ${p1_files} | sort -V)
    p2_sorted=\$(printf '%s\\n' ${p2_files} | sort -V)

    cat \$p1_sorted > all_R1.fastq.gz
    cat \$p2_sorted > all_R2.fastq.gz

    for f in \$p1_sorted; do
        sample_id=\$(basename \$f | cut -d'.' -f1)
    done
    paste <(printf '%s\\n' \$p1_sorted) <(printf '%s\\n' \$p2_sorted) | awk -F'\\t' '{split(\$1,a,"/"); split(a[length(a)],b,"."); print b[1]"\t"\$1"\t"\$2}' > manifest.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        script: "merge_corrected.nf"
    END_VERSIONS
    """
}
