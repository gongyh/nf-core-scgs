process MARKER_NCLUSTERS {
    tag "${fasta_file.baseName}"

    conda "${moduleDir}/copygen.yaml"
    container 'community.wave.seqera.io/library/biopython_einops_loguru_numpy_pruned:f627fd125f815a95'

    input:
    path kmer_file
    path fasta_file

    output:
    path "cluster_value", emit: marker_cv

    script:
    def args    = task.ext.args ?: ''
    """
    python -c "import sys; print('Python version:', sys.version); print('sys.path:', sys.path)"
    python -c "import sklearn; print('sklearn version:', sklearn.__version__)"
    python ${projectDir}/bin/dcvbin/marker_gene/src/marker_gene_utils.py \
        -kf "${kmer_file}" \
        -cf "${fasta_file}" \
        -cvf "cluster_value"
    """
}
