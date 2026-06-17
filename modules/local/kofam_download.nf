process KOFAM_DOWNLOAD {
    tag "${db_url}"
    publishDir "${params.outdir}/kofam", mode: 'copy'

    input:
    val(db_url)
    val(out_dir)

    output:
    path 'kofam_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p kofam_db
    cd kofam_db
    
    echo "Downloading KOfam database from ${db_url}..."
    
    wget -q "${db_url}profiles.tar.gz" && tar -xzf profiles.tar.gz && rm profiles.tar.gz
    wget -q "${db_url}ko_list"
    
    echo "KOfam database downloaded successfully"
    
    cat > versions.yml << 'EOF'
"kofam_download":
  "version": "1.0.0"
EOF
    """
}
