process MARKER_NCLUSTERS {
    tag "${meta.id}"

    conda "${moduleDir}/copygen.yaml"
    container 'community.wave.seqera.io/library/copygen:aca96b4a00a56131'

    input:
    path kmer_file
    path fasta_file

    output:
    path "cluster_value", emit: marker_cv

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/dcvbin/marker_gene/src/marker_gene_utils.py \
        -kf "${kmer_file}" \
        -cf "${fasta_file}" \
        -cvf "cluster_value"
    """
}
