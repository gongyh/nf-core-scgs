process BOWTIE2_ALIGN {
    tag "$prefix"
    label "process_high"
    publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/bowtie2" : params.outdir }, mode: 'copy',
                saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }

    conda "bioconda::bowtie2=2.4.4 bioconda::samtools=1.16.1 conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0' :
        'biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0' }"

    input:
    path reads
    path index
    val(single_end)

    output:
    path '*.bam', emit: bb_bam

    when:
    denovo == false

    script:
    prefix = reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    filtering = params.allow_multi_align ? '' : "| samtools view -b -q 40 -F 4 -F 256 -"
    R1 = reads[0].toString()
    if (single_end) {
    """
    bowtie2 -x ${index}/genome -p ${task.cpus} -U $R1 | samtools view -bT $index - $filtering > ${prefix}.bam
    """
    } else {
    R2 = reads[1].toString()
    """
    bowtie2 --no-mixed --no-discordant -X 1000 -x ${index}/genome -p ${task.cpus} -1 $R1 -2 $R2 | samtools view -bT $index - $filtering > ${prefix}.bam
    """
    }
}
