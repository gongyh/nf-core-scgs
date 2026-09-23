nextflow.enable.types = true

process CONTIG_COVERAGE {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple(meta: Map, bam: Path, bai: List<Path>, fasta: Path, fai: Path)

    output:
    record(meta: meta, depth: file("${meta.id}.depth"), mqc_tsv: file("${meta.id}_coverage_mqc.tsv"))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    if [ ! -f "${bam}.bai" ]; then
        samtools index "${bam}"
    fi
    samtools coverage --reference "${fasta}" -o "${meta.id}.cov" "${bam}"
    awk '!/^#/ {print \$1"\\t"\$7}' "${meta.id}.cov" | sort -k1,1 > "${meta.id}.depth"
    rm "${meta.id}.cov"

    cat > "${meta.id}_coverage_mqc.tsv" <<'EOF'
# id: coverage
# section_name: Coverage
# plot_type: table
EOF
    printf 'Sample\\tContigs\\tContigs with coverage\\tMean depth\\n' >> "${meta.id}_coverage_mqc.tsv"
    awk -F '\\t' -v sample='${meta.id}' '
        { contigs++; if (\$2 > 0) covered++; depth += \$2 }
        END { printf "%s\\t%d\\t%d\\t%.4f\\n", sample, contigs, covered, contigs ? depth / contigs : 0 }
    ' "${meta.id}.depth" >> "${meta.id}_coverage_mqc.tsv"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/^.*samtools //')
    END_VERSIONS
    """
}

process MERGE_COVERAGE {
    tag "merge_all"
    label 'process_low'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    depth_files: Bag<Path>

    output:
    record(matrix: file('abundance_matrix.tsv'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    samples=(\$(ls *.depth | sed 's/.depth//'))
    cut -f1 *.depth | sort -u > all_contigs.tmp
    for sample in \${samples[*]}; do
        join -a1 -e0 -o '2.2' -t \$'\\t' all_contigs.tmp "\${sample}.depth" > "\${sample}.depth_col"
    done
    paste all_contigs.tmp \$(for s in \${samples[*]}; do echo "\${s}.depth_col"; done) > abundance_matrix.tsv
    header="contig_id"
    for sample in \${samples[*]}; do
        header="\${header}\\t\${sample}"
    done
    (echo -e "\${header}" && cat abundance_matrix.tsv) > abundance_matrix.tsv.tmp && mv abundance_matrix.tsv.tmp abundance_matrix.tsv
    rm -f *.depth *.depth_col all_contigs.tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge: bash \$(bash --version | head -1)
    END_VERSIONS
    """
}
