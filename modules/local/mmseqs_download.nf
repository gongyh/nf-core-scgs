process MMSEQS_DOWNLOAD {
    tag "${db_url}"
    publishDir "${params.outdir}/mmseqs", mode: 'copy'

    input:
    val(db_url)
    val(out_dir)

    output:
    path 'mmseqs_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p mmseqs_db
    cd mmseqs_db
    
    echo "Downloading MMseqs2 database from ${db_url}..."
    
    wget -q -r -np -nH --cut-dirs=1 -R "index.html*" "${db_url}" || \\\n    curl -L -o mmseqs_db.tar.gz "${db_url}" && tar -xzf mmseqs_db.tar.gz
    
    echo "MMseqs2 database downloaded successfully"
    
    cat > versions.yml << 'EOF'
"mmseqs_download":
  "version": "1.0.0"
EOF
    """
}
