process CHECKM2_DBDOWNLOAD {
    tag "CheckM2"

    conda "bioconda::checkm2=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0' :
        'community.wave.seqera.io/library/checkm2:1.0.1--034a3a15afae63b1' }"

    output:
    path 'checkm2_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p checkm2_db
    echo "Downloading CheckM2 database ..."
    checkm2 database --download --path checkm2_db
    echo "CheckM2 database downloaded successfully"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$( checkm2 --version )
    END_VERSIONS
    """
}
