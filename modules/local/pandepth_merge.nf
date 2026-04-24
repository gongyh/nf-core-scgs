process PANDEPTH_MERGE {
    tag "merge_all"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    path depth_files
    output:
    path "abundance_matrix.tsv", emit: matrix
    path "versions.yml"        , emit: versions

    script:
    """
    samples=(\$(ls *.depth | sed 's/.depth//'))
    cut -f1 *.depth | sort -u > all_contigs.tmp
    for sample in \${samples[*]}; do
        join -a1 -e0 -o '2.2' -t \$'\t' all_contigs.tmp "\${sample}.depth" > "\${sample}.depth_col"
    done
    paste all_contigs.tmp \$(for s in \${samples[*]}; do echo "\${s}.depth_col"; done) > abundance_matrix.tsv
    header="contig_id"
    for sample in \${samples[*]}; do
        header="\${header}\t\${sample}"
    done
    (echo -e "\${header}" && cat abundance_matrix.tsv) > abundance_matrix.tsv.tmp && mv abundance_matrix.tsv.tmp abundance_matrix.tsv
    rm -f *.depth *.depth_col all_contigs.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge: bash \$(bash --version | head -1)
    END_VERSIONS
    """
}
