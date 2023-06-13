process RESFINDER {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(contigs)
    path db

    output:
    path("${prefix}/*")

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $prefix
    python /opt/resfinder/resfinder.py -i $contigs -o $prefix -p $db -mp blastn -x
    rm -rf $prefix/tmp
    """
}
