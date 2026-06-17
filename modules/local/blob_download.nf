process BLOB_DOWNLOAD {
    tag "${db_url}"
    publishDir "${params.outdir}/blob", mode: 'copy'

    input:
    val(db_url)
    val(out_dir)

    output:
    path 'blob_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p blob_db
    cd blob_db

    echo "Downloading Blobtools nodesDB from ${db_url}..."

    wget -q "${db_url}/taxdump.tar.gz" && tar -xzf taxdump.tar.gz && rm taxdump.tar.gz

    echo "Blobtools database downloaded successfully"

    cat > versions.yml << 'EOF'
"blob_download":
    "version": "1.0.0"
EOF
    """
}
