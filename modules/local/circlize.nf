process CIRCLIZE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.0=h468198e_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3' :
        'biocontainers/bedtools:2.31.0--h468198e_0' }"

    input:
    tuple val(meta), path(sbed)
    path(refbed)

    output:
    tuple val(meta), path("${prefix}-cov200.bed"), emit: bed
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools makewindows -b $refbed -w 200 > genome.200.bed
    bedtools coverage -mean -b $sbed -a genome.200.bed | sort -k 1,1 -V -k 2n,2 -k 3n,3 > ${prefix}-cov200.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//; s/Using.*\$//')
    END_VERSIONS
    """
}
