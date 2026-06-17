process KRAKEN2_DOWNLOAD {
    tag "${db_url}"
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    val(db_url)
    val(out_dir)

    output:
    path 'kraken2_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p kraken2_db
    cd kraken2_db

    echo "Downloading Kraken2 database from ${db_url}..."

    wget -q "${db_url}/k2_standard_08gb_20231009.tar.gz" && \\\n    tar -xzf k2_standard_08gb_20231009.tar.gz && rm k2_standard_08gb_20231009.tar.gz

    echo "Kraken2 database downloaded successfully"

    cat > versions.yml << 'EOF'
"kraken2_download":
    "version": "1.0.0"
EOF
    """
}
