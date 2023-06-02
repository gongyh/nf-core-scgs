process SAMTOOLS {
    tag "${prefix}"
    label 'process_medium'
    publishDir path: "${pp_outdir}", mode: 'copy',
                saveAs: { filename ->
                    if (filename.indexOf(".stats.txt") > 0) "stats/$filename"
                    else if (filename.indexOf("_bins.txt") > 0) filename
                    else if (filename.indexOf("_pdrc.pdf") > 0) filename
                    else params.saveAlignedIntermediates ? filename : null
                }

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    path bam
    path genome

    output:
    path '*.markdup.bam',                      emit: bam
    path '*.markdup.bam.bai',                  emit: bai
    path '*.markdup.bed',                      emit: bed
    path '*.stats.txt',                        emit: stats
    path "${prefix}_1k_bins.txt"
    path "${prefix}_pdrc.pdf"

    script:
    pp_outdir = "${params.outdir}/bowtie2"
    prefix = bam.baseName
    """
    samtools sort -o ${prefix}.sorted.bam $bam
    samtools index ${prefix}.sorted.bam
    picard MarkDuplicates I=${prefix}.sorted.bam O=${prefix}.markdup.bam M=metrics.txt AS=true
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
    """
}
