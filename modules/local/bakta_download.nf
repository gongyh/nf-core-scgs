nextflow.enable.types = true

process BAKTA_DBDOWNLOAD {
    tag 'Bakta'
    label 'process_single'

    conda "bioconda::bakta=1.11.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.11.4--pyhdfd78af_0' :
        'biocontainers/bakta:1.11.4--pyhdfd78af_0' }"

    output:
    record(db: file('bakta_db', type: 'dir'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    bakta_db download --output bakta_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(bakta_db --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bakta_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: '1.11.4'
    END_VERSIONS
    """
}
