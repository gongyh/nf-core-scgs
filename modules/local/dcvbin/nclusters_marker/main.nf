process MARKER_NCLUSTERS {
    tag "${fasta_file.baseName}"

    conda "${moduleDir}/copygen.yaml"
    container 'community.wave.seqera.io/library/numpy_pandas_python-abi3:f58d0c4ace38e4d4'

    input:
    path kmer_file
    path fasta_file

    output:
    path "cluster_value", emit: marker_cv

    script:
    def args    = task.ext.args ?: ''
    """
    python ${projectDir}/bin/dcvbin/marker_gene/src/marker_gene_utils.py \
        -kf "${kmer_file}" \
        -cf "${fasta_file}" \
        -cvf "cluster_value"
    """
}
