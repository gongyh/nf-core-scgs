nextflow.enable.types = true

process GTDB_DBDOWNLOAD {
    tag "GTDB_r214"

    conda "conda-forge::wget=1.25.0"
    container "community.wave.seqera.io/library/wget:1.25.0--817c089a96769e94"

    output:
    record(db: file('gtdb_db'), versions: file('versions.yml'))

    script:
    """
    mkdir -p gtdb_db
    echo "Downloading GTDB database ..."
    wget -q https://data.gtdb.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
    tar xvzf gtdbtk_r232_data.tar.gz -C gtdb_db/
    rm -f gtdbtk_r232_data.tar
    echo "GTDB database downloaded successfully"

    cat > versions.yml << 'EOF'
    "gtdb_download":
        "version": "1.0.0"
    EOF
    """
}
