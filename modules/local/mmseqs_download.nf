process MMSEQS_DBDOWNLOAD {
    tag "${MMseqs2}"

    conda "bioconda::mmseqs2=18.8cc5c conda-forge::wget=1.25.0"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/edfecaaca16ca7fb7b6428dce0ed9c737549b38146360c98fdabf74e6c4cac68/data'
        : 'community.wave.seqera.io/library/mmseqs2_wget:aa683a2c5355899d'}"

    output:
    path 'mmseqs_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p mmseqs_db tmp
    echo "Downloading MMseqs2 database ..."
    mmseqs databases GTDB mmseqs_db tmp --threads ${task.cpus}
    rm -rf tmp
    echo "MMseqs2 database downloaded successfully"

    cat > versions.yml << 'EOF'
"mmseqs_download":
    "version": "1.0.0"
EOF
    """
}
