nextflow.enable.types = true

process METABULI_DBDOWNLOAD {
    tag "Metabuli_GTDB226"

    conda "conda-forge::wget=1.25.0"
    container "community.wave.seqera.io/library/wget:1.25.0--817c089a96769e94"

    output:
    record(db: file('metabuli_db'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    mkdir -p metabuli_db
    echo "Downloading MetaBuli database ..."
    wget -q "https://opendata.mmseqs.org/metabuli/gtdb226.tar.gz"
    tar -xzf gtdb226.tar.gz -C metabuli_db/ && rm -f gtdb226.tar.gz
    echo "MetaBuli database downloaded successfully"

    cat > versions.yml << 'EOF'
"metabuli_download":
    "version": "1.0.0"
EOF
    """

    stub:
    """
    mkdir -p metabuli_db

    cat > versions.yml << 'EOF'
"metabuli_download":
    "version": "1.0.0"
EOF
    """
}
