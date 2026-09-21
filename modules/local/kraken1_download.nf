nextflow.enable.types = true

process KRAKEN1_DBDOWNLOAD {
    tag 'Kraken1 MiniKraken'
    label 'process_single'

    conda "conda-forge::wget=1.25.0"
    container "community.wave.seqera.io/library/wget:1.25.0--817c089a96769e94"

    output:
    record(db: file('kraken1_db', type: 'dir'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    wget -q https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz
    tar -xzf minikraken_20171019_8GB.tgz
    mv minikraken_20171019_8GB kraken1_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken1_database: 'minikraken_20171019_8GB'
    END_VERSIONS
    """

    stub:
    """
    mkdir -p kraken1_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken1_database: 'minikraken_20171019_8GB'
    END_VERSIONS
    """
}
