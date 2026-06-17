process EGGNOG_DBDOWNLOAD {
    tag "eggNOG"

    conda "bioconda::eggnog-mapper=2.1.11=pyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.11--pyhdfd78af_0' :
        'biocontainers/eggnog-mapper:2.1.11--pyhdfd78af_0' }"

    output:
    path 'eggnog_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p eggnog_db
    echo "Downloading EggNOG database ..."
    download_eggnog_data.py --data_dir eggnog_db
    echo "EggNOG database downloaded successfully"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$(echo \$(emapper.py --version | grep emapper 2>&1 ) | cut -d'/' -f1 | sed 's/^.*emapper-//; s/Using.*\$//')
    END_VERSIONS
    """
}
