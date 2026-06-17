process GTDB_DOWNLOAD {
    tag "${db_url}"
    publishDir "${params.outdir}/gtdb", mode: 'copy'

    input:
    val(db_url)
    val(out_dir)

    output:
    path 'gtdb_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p gtdb_db
    cd gtdb_db

    echo "Downloading GTDB database from ${db_url}..."

    wget -q -r -np -nH --cut-dirs=3 -R "index.html*" "${db_url}/latest/" || \\\n    wget -q "${db_url}/release214/auxillary_files/gtdbtk_r214_data.tar.gz" && tar -xzf gtdbtk_r214_data.tar.gz

    echo "GTDB database downloaded successfully"

    cat > versions.yml << 'EOF'
"gtdb_download":
    "version": "1.0.0"
EOF
    """
}
