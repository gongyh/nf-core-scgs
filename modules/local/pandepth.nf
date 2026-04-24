process CONTIG_COVERAGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.depth"), emit: depth
    path "versions.yml", emit: versions

    script:
    def pandepth_bin = "${projectDir}/bin/pandepth"
    def threads = task.cpus ?: 1
    """
    if [ ! -f "${bam}.bai" ]; then
        samtools index "${bam}"
    fi
    ${pandepth_bin} -i "${bam}" -t ${threads} -r "${fasta}" -o "${meta.id}"
    zcat "${meta.id}.chr.stat.gz" 2>/dev/null | awk 'NR>1 {print \$1"\\t"\$5}' | sort -k1,1 > "${meta.id}.depth"
    rm -f "${meta.id}.chr.stat.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandepth: \$(${pandepth_bin} -h 2>&1 | head -1)
    END_VERSIONS
    """
}

process MERGE_COVERAGE {
    tag "merge_all"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5"

    input:
    path depth_files

    output:
    path "abundance_matrix.tsv", emit: matrix
    path "versions.yml", emit: versions

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
