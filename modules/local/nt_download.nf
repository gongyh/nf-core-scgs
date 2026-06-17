process NT_DOWNLOAD {
    tag "${db_url}"
    publishDir "${params.outdir}/nt", mode: 'copy'

    input:
    val(db_url)
    val(out_dir)

    output:
    path 'nt_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p nt_db
    cd nt_db

    echo "Downloading NCBI nt database from ${db_url}..."

    wget -q "${db_url}/nt.00.tar.gz" && tar -xzf nt.00.tar.gz && rm nt.00.tar.gz

    echo "NCBI nt database downloaded successfully"

    cat > versions.yml << 'EOF'
"nt_download":
    "version": "1.0.0"
EOF
    """
}
