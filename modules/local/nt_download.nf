process NT_DBDOWNLOAD {
    tag "NCBI_nt"

    conda "bioconda::blast=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.13.0--hf3cf87c_0' :
        'biocontainers/blast:2.13.0--hf3cf87c_0' }"

    output:
    path 'nt_db', emit: db
    path 'versions.yml', emit: versions

    script:
    """
    mkdir -p nt_db
    cd nt_db
    echo "Downloading NCBI nt database ..."
    update_blastdb.pl --decompress nt
    echo "NCBI nt database downloaded successfully"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(blastn -version 2>&1 | grep blastn | sed 's/^.*blastn: //; s/Using.*\$//')
    END_VERSIONS
    """
}
