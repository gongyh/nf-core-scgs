process METABULI_DOWNLOAD {
    tag "${db_url}"
    publishDir "${params.outdir}/metabuli", mode: 'copy'

    input:
    val(db_url)
    val(out_dir)

    output:
    path 'metabuli_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p metabuli_db
    cd metabuli_db

    echo "Downloading MetaBuli database from ${db_url}..."

    git clone --depth 1 "${db_url}" . || \\\n    wget -q "https://github.com/khyox/metabuli/releases/download/v1.0/metabuli_db.tar.gz" && tar -xzf metabuli_db.tar.gz

    echo "MetaBuli database downloaded successfully"

    cat > versions.yml << 'EOF'
"metabuli_download":
    "version": "1.0.0"
EOF
    """
}
