process KOFAM_DBDOWNLOAD {
    tag "${kofam}"

    conda "conda-forge::wget=1.25.0"
    container "community.wave.seqera.io/library/wget:1.25.0--817c089a96769e94"

    output:
    path 'kofam_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p kofam_db
    echo "Downloading KOfam database from https://www.genome.jp/ftp/db/kofam/ ..."
    wget -q "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz"
    tar -xzf profiles.tar.gz -C kofam_db/ && rm -f profiles.tar.gz
    wget -q "https://www.genome.jp/ftp/db/kofam/ko_list.gz"
    gzip -cd ko_list.gz > kofam_db/ko_list && rm -f ko_list.gz
    echo "KOfam database downloaded successfully"

    cat > versions.yml << 'EOF'
"kofam_download":
    "version": "1.0.0"
EOF
    """
}
