process QUALIMAP {
    label 'process_medium'
    publishDir "${pp_outdir}", mode: 'copy'

    conda "bioconda::qualimap=2.2.2d"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--1' :
        'biocontainers/qualimap:2.2.2d--1' }"

    input:
    path("*")
    path("*")
    path gff

    output:
    path('*.markdup_stats'),                  emit: qualimap
    path('multi-bamqc')

    script:
    pp_outdir = "${params.outdir}/qualimap"
    """
    ls *.markdup.bam > bams.txt
    let num=`ls *.bam | wc -l`
    if [ \$num == 1 ]; then
    qualimap bamqc -c -bam *.markdup.bam -gff $gff -outdir multi-bamqc
    ln -s multi-bamqc Sample.markdup_stats
    else
    cat bams.txt | awk '{split(\$1,a,".markdup.bam"); print a[1]"\t"\$1}' > inputs.txt
    JAVA_MEM_SIZE=${task.memory.toGiga()}G qualimap multi-bamqc -r -c -d inputs.txt -gff $gff -outdir multi-bamqc
    fi
    """
}
