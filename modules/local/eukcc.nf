process EUKCC {
    publishDir "${params.outdir}/EukCC", mode: 'copy'

    input:
    path faa
    path db

    output:
    path("${prefix}")

    when:
    euk && eukcc_db

    script:
    prefix = faa.toString() - ~/(\.faa)?(\.aa)?(\.fasta)?(\.fa)?$/
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
