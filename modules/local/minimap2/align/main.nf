process MINIMAP2_ALIGN {
    tag "$prefix"
    label 'process_medium'
    publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/bowtie2" : params.outdir }, mode: 'copy',
                saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "bioconda::minimap2=2.24 bioconda::samtools=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    path reads
    path fasta

    output:
    path '*.bam', emit: bb_bam

    when:
    denovo == false

    script:
    prefix = reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    filtering = params.allow_multi_align ? '' : "-F 256"
    R1 = reads[0].toString()
    """
    minimap2 -x map-ont -a $fasta $R1 | samtools view -b -q 40 -F 4 $filtering - > ${prefix}.bam
    """
}
