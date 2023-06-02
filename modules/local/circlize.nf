process CIRCLIZE {
    tag "${prefix}"
    label 'process_medium'
    publishDir "${params.outdir}/circlize", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bed") > 0) "$filename" else null
            }

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3':
        'biocontainers/bedtools:2.30.0--h468198e_3' }"

    input:
    path sbed
    path refbed

    output:
    path("${prefix}-cov200.bed")

    script:
    prefix = sbed.toString() - ~/(\.markdup\.bed)?(\.markdup)?(\.bed)?$/
    """
    bedtools makewindows -b $refbed -w 200 > genome.200.bed
    bedtools coverage -mean -b $sbed -a genome.200.bed | sort -k 1V,1 -k 2n,2 -k 3n,3 > ${prefix}-cov200.bed
    """
}
