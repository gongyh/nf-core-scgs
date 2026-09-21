nextflow.enable.types = true

process BLOB_DBDOWNLOAD {
    tag "BlobTools"

    conda "bioconda::blobtools=1.1.1=py_1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blobtools:1.1.1--py_1' :
        'biocontainers/blobtools:1.1.1--py_1' }"

    output:
    record(db: file('blob_db'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    mkdir -p blob_db data
    echo "Downloading Blobtools nodesDB ..."
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
    tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
    blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp
    cp \$(echo \$(which blobtools))/data/nodesDB.txt blob_db/nodesDB.txt
    echo "Blobtools database downloaded successfully"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools: \$(echo \$(blobtools -v 2>&1) | sed 's/^.*blobtools v//; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p blob_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools: '1.1.1'
    END_VERSIONS
    """
}
