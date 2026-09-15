process GENOMAD_DBDOWNLOAD {
    tag "geNomad"

    conda "bioconda::genomad=1.7.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.7.4--pyhdfd78af_0':
        'community.wave.seqera.io/library/genomad:1.7.4--605ab516f999b1b4' }"

    output:
    path 'db'          , emit: db
    path 'versions.yml', emit: versions

    script:
    """
    echo "Downloading GENOMAD database ..."
    genomad download-database .
    cp -r genomad_db db
    rm -rf genomad_db
    echo "GENOMAD database downloaded successfully"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """
}
