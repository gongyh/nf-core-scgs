process SAMTOOLS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::samtools=1.17 bioconda::picard=2.19.0 bioconda::bedtools=2.31.0 conda-forge::r-magicaxis=2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a055626708da8b97b79d76d03667dba1cb4bb107:3785590b864aff1747a8d3237666cee5030d4c3e-0' :
        'scgs/mulled-v2-a055626708da8b97b79d76d03667dba1cb4bb107:3785590b864aff1747a8d3237666cee5030d4c3e-0' }"

    input:
    tuple val(meta), path(bam)
    path genome

    output:
    tuple val(meta), path("*.markdup.bam")    , emit: bam
    tuple val(meta), path("*.markdup.bam.bai"), emit: bai
    tuple val(meta), path("*.markdup.bed")    , emit: bed
    tuple val(meta), path("*.stats.txt")      , optional:true, emit: stats
    path  "versions.yml"                      , emit: versions
    path("${prefix}_1k_bins.txt")             , emit: txt
    path("${prefix}_pdrc.pdf")                , emit: pdf

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort -o ${prefix}.sorted.bam $bam
    samtools index ${prefix}.sorted.bam
    picard -Xmx${task.memory.toGiga()}G MarkDuplicates I=${prefix}.sorted.bam O=${prefix}.markdup.bam M=metrics.txt AS=true
    samtools index ${prefix}.markdup.bam
    bedtools bamtobed -i ${prefix}.markdup.bam | sort -T /tmp -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${prefix}.markdup.bed
    samtools stats -t ${genome} ${prefix}.markdup.bam > ${prefix}.stats.txt
    # uniformity
    cut -f1,3 ${genome} > ref.genome
    genomeCoverageBed -ibam ${prefix}.markdup.bam -d -g ref.genome > ${prefix}.raw.cov
    meanCov=\$(awk 'BEGIN{ total=0; base=0 } { total=total+\$3; base=base+1 } END{ printf total/base }' ${prefix}.raw.cov)
    awk -v mc=\$meanCov -F'\t' '{print \$1"\t"\$2"\t"\$3/mc}' ${prefix}.raw.cov > ${prefix}.relative.cov
    awk '{sum+=\$3} (NR%1000)==0{print sum/1000; sum=0;}' ${prefix}.relative.cov > ${prefix}_1k_bins.txt
    plotProp.R ${prefix}_1k_bins.txt ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//; s/Using.*\$//')
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1 | grep SNAPSHOT | cut -d'-' -f1))
    END_VERSIONS
    """
}
