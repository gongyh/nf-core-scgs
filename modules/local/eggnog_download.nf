process EGGNOG_DOWNLOAD {
    tag "${db_url}"
    publishDir "${params.outdir}/eggnog", mode: 'copy'

    input:
    val(db_url)
    val(out_dir)

    output:
    path 'eggnog_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p eggnog_db
    cd eggnog_db

    echo "Downloading EggNOG database from ${db_url}..."

    wget -q "${db_url}eggnog.db.gz" && gunzip eggnog.db.gz
    wget -q "${db_url}eggnog.taxid_info.tsv.gz" && gunzip eggnog.taxid_info.tsv.gz
    wget -q "${db_url}members.tsv.gz" && gunzip members.tsv.gz

    echo "EggNOG database downloaded successfully"

    cat > versions.yml << 'EOF'
"eggnog_download":
    "version": "1.0.0"
EOF
    """
}
