process ACDC {
    tag "$meta.id"

    input:
    tuple val(meta), path contigs
    tuple val(meta), path tax
    path db

    output:
    path("${prefix}")

    when:
    false

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $tax | grep -v '^#' | cut -f1,18 > genus.txt
    /usr/local/bin/acdc -i $contigs -m 1000 -b 100 -o $prefix -K $db -x genus.txt -T ${task.cpus}
    """
}
