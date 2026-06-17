process GENOMAD_DOWNLOAD {
    tag "${db_url}"
    publishDir "${params.outdir}/genomad", mode: 'copy'

    input:
    val(db_url)
    val(out_dir)

    output:
    path 'genomad_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p genomad_db
    cd genomad_db

    echo "Downloading GENOMAD database from ${db_url}..."

    wget -q "${db_url}/viral_db.tar.gz" && tar -xzf viral_db.tar.gz && rm viral_db.tar.gz

    echo "GENOMAD database downloaded successfully"

    cat > versions.yml << 'EOF'
"genomad_download":
    "version": "1.0.0"
EOF
    """
}
