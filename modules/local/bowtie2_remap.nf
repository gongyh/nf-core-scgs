process BOWTIE2_REMAP {
    tag "$meta.id"

    input:
    tuple val(meta), path(contigs)

    output:
    path("${prefix}Bowtie2Index"),    emit: index

    when:
    params.remap

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}Bowtie2Index; cd ${prefix}Bowtie2Index
    ln -s ../${contigs} ${prefix}.fa
    bowtie2-build ${prefix}.fa ${prefix}
    """
}
