process EUKCC {
    tag "$meta.id"
    label 'process_medium'

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
