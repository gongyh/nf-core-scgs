process DCVBIN_BIN {
    tag "${meta.id}"
    
    conda "${moduleDir}/dnaberts.yaml"
    container 'community.wave.seqera.io/library/dnaberts:7a7299083f265248'

    input:
    path vae_features_file
    path cluster_value_file
    path fasta_file
    
    output:
    path "${prefix}_bins",        emit: bins_dir
    path "${prefix}_prinum.txt",  emit: label_file
    
    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    export QT_XCB_GL_INTEGRATION="none"
    mkdir -p "${prefix}_bins"
    python ${projectDir}/bin/dcvbin/scripts/myCluster_2.py \
        -vd "${vae_features_file}" \
        -ld "${prefix}_prinum.txt" \
        -fd "${fasta_file}" \
        -bd "${prefix}_bins" \
        -cvf "${cluster_value_file}"
    """
}
