process KRAKEN2_DOWNLOAD {
    tag "${Kraken2}"

    conda "conda-forge::wget=1.25.0"
    container "community.wave.seqera.io/library/wget:1.25.0--817c089a96769e94"

    output:
    path 'kraken2_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p kraken2_db
    echo "Downloading Kraken2 database ..."
    wget -q "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260226.tar.gz"
    tar -xzf k2_pluspf_20260226.tar.gz -C kraken2_db/ && rm -f k2_pluspf_20260226.tar.gz
    echo "Kraken2 database downloaded successfully"

    cat > versions.yml << 'EOF'
"kraken2_download":
    "version": "1.0.0"
EOF
    """
}
