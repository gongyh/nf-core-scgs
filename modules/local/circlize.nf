process CIRCLIZE {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/circlize", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bed") > 0) "$filename" else null
            }

    input:
    tuple val(meta), path(sbed)
    path refbed

    output:
    path("${prefix}-cov200.bed")

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bedtools makewindows -b $refbed -w 200 > genome.200.bed
    bedtools coverage -mean -b $sbed -a genome.200.bed | sort -k 1V,1 -k 2n,2 -k 3n,3 > ${prefix}-cov200.bed
    """
}
