nextflow.enable.types = true

process DNABERTS_DBDOWNLOAD {
    tag 'DNABERT-S'
    label 'process_single'

    conda "conda-forge::wget=1.25.0"
    container "community.wave.seqera.io/library/wget:1.25.0--817c089a96769e94"

    output:
    record(db: file('dnaberts_db', type: 'dir'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    wget -q -O DNABERT-S.tar.gz \
        https://github.com/gongyh/dcvbin/releases/download/v1.1/DNABERT-S.tar.gz
    mkdir -p dnaberts_db
    tar -xzf DNABERT-S.tar.gz --strip-components=1 -C dnaberts_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnaberts_model: 'Figshare article 30400258'
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dnaberts_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnaberts_model: 'Figshare article 30400258'
    END_VERSIONS
    """
}
