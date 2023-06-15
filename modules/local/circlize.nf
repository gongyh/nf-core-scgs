process CIRCLIZE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(sbed)
    path(refbed)

    output:
    path("${prefix}-cov200.bed")
    path "versions.yml",  emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools makewindows -b $refbed -w 200 > genome.200.bed
    bedtools coverage -mean -b $sbed -a genome.200.bed | sort -k 1V,1 -k 2n,2 -k 3n,3 > ${prefix}-cov200.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
