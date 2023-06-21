process EUKCC {
    tag "$meta.id"
    label 'process_medium'

    conda "eukcc=2.1.0--pypyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukcc:2.1.0--pypyhdfd78af_0' :
        'biocontainers/eukcc:2.1.0--pypyhdfd78af_0' }"

    input:
    tuple val(meta), path(faa)
    path db

    output:
    path("${prefix}")

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export HOME=/tmp/
    if [ -f "/tmp/.etetoolkit/taxa.sqlite" ]; then
    echo "NCBI taxa database exist!"
    else
    python -c "from ete3 import NCBITaxa; ncbi = NCBITaxa(taxdump_file='/opt/nf-core-scgs/taxdump.tar.gz')"
    fi
    eukcc --db ${db} --ncores ${task.cpus} --outdir ${prefix} --protein ${faa} || echo "Ignore minor errors of eukcc!"
    """
}
