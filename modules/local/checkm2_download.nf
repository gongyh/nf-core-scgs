process CHECKM2_DOWNLOAD {
    tag "${db_url}"
    publishDir "${params.outdir}/checkm2", mode: 'copy'

    input:
    val(db_url)
    val(out_dir)

    output:
    path 'checkm2_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p checkm2_db
    cd checkm2_db

    echo "Downloading CheckM2 database from ${db_url}..."

    wget -q -r -np -nH --cut-dirs=2 -R "index.html*" "${db_url}" || \\\n    curl -L -o checkm2_db.tar.gz "${db_url}checkm2_database.tar.gz" && tar -xzf checkm2_db.tar.gz

    echo "CheckM2 database downloaded successfully"

    cat > versions.yml << 'EOF'
"checkm2_download":
  "version": "1.0.0"
EOF
    """
}
