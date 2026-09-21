nextflow.enable.types = true

process EUKCC_DBDOWNLOAD {
    tag 'EukCC'
    label 'process_single'

    conda "bioconda::eukcc=2.1.0=pypyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukcc:2.1.0--pypyhdfd78af_0' :
        'biocontainers/eukcc:2.1.0--pypyhdfd78af_0' }"

    output:
    record(db: file('eukcc_db', type: 'dir'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    python - <<'PY'
    from urllib.request import urlretrieve

    urlretrieve(
        'http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.2.tar.gz',
        'eukcc2_db_ver_1.2.tar.gz',
    )
    PY
    tar -xzf eukcc2_db_ver_1.2.tar.gz
    mv eukcc2_db_ver_1.2 eukcc_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eukcc: \$(eukcc -v 2>&1 | sed 's/^.*EukCC version //; s/Using.*\$//')
        eukcc_database: '1.2'
    END_VERSIONS
    """

    stub:
    """
    mkdir -p eukcc_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eukcc: '2.1.0'
        eukcc_database: '1.2'
    END_VERSIONS
    """
}
