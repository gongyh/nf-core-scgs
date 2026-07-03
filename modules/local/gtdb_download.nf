process GTDB_DBDOWNLOAD {
    tag "GTDB_r214"

    conda "conda-forge::wget=1.25.0"
    container "community.wave.seqera.io/library/wget:1.25.0--817c089a96769e94"

    output:
    path 'gtdb_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p gtdb_db
    echo "Downloading GTDB database ..."
    wget -q https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz
    tar xvzf gtdbtk_r214_data.tar.gz -C gtdb_db/
    rm -f gtdbtk_r214_data.tar.gz
    echo "GTDB database downloaded successfully"

    cat > versions.yml << 'EOF'
    "gtdb_download":
        "version": "1.0.0"
    EOF
    """
}
