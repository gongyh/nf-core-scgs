process MERGE_CORRECTED {
    tag "merge"
    label 'process_low'
    publishDir "${params.outdir}/merged", mode: 'copy'
    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"
    input:
    path p1_files
    path p2_files

    output:
    path "all_R1.fastq.gz", emit: r1
    path "all_R2.fastq.gz", emit: r2
    path "manifest.txt", emit: manifest
    path "manifest_mqc.tsv", emit: manifest_mqc

    script:
    """
    p1_sorted=\$(printf '%s\\n' ${p1_files} | sort -V)
    p2_sorted=\$(printf '%s\\n' ${p2_files} | sort -V)
    cat \$p1_sorted > all_R1.fastq.gz
    cat \$p2_sorted > all_R2.fastq.gz

    cat > manifest_mqc.tsv <<'EOF'
# id: "merged_samples_metadata"
# section_name: "Sample Metadata"
# description: "List of merged samples with their original read file paths."
# format: "tsv"
# plot_type: "table"
# Sample_ID	R1_File_Path	R2_File_Path
EOF

    > manifest.txt

    paste <(printf '%s\\n' \$p1_sorted) <(printf '%s\\n' \$p2_sorted) | while IFS=\$'\t' read r1 r2;
do
        id=\$(basename "\$r1" | cut -d'.' -f1)
        printf "%s\\t%s\\t%s\\n" "\$id" "\$r1" "\$r2" >> manifest.txt
        printf "%s\\t%s\\t%s\\n" "\$id" "\$r1" "\$r2" >> manifest_mqc.tsv
    done
    """
}
