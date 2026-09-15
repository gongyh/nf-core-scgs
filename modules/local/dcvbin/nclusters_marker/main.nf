process MARKER_NCLUSTERS {
    tag "$meta.id"

    conda "${moduleDir}/copygen.yaml"
    container 'community.wave.seqera.io/library/copygen:aca96b4a00a56131'

    input:
    tuple val(meta), path(kmer_file)
    tuple val(meta), path(fasta_file)

    output:
    tuple val(meta), path("cluster_value"), emit: marker_cv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/dcvbin/marker_gene/src/marker_gene_utils.py \
        -kf "${kmer_file}" \
        -cf "${fasta_file}" \
        -cvf "cluster_value"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copygen: \$(python -c "import copygen; print(copygen.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
