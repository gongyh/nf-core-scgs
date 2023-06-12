process INDELREALIGN {
    tag "${meta.id}"
    label 'process_single'
    publishDir "${pp_outdir}", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path fa

    output:
    tuple val(meta), path("*.realign.bam"),         emit: bam
    tuple val(meta), path("*.realign.bam.bai"),     emit: bai

    when:
    params.snv && !params.nanopore

    script:
    pp_outdir = "${params.outdir}/gatk"
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    samtools faidx $fa
    picard CreateSequenceDictionary R=$fa
    picard AddOrReplaceReadGroups I=$bam O=${prefix}_rg.bam RGLB=lib RGPL=illumina RGPU=run RGSM=${prefix}
    samtools index ${prefix}_rg.bam
    gatk3 -T RealignerTargetCreator -R $fa -I ${prefix}_rg.bam -o indels.intervals
    gatk3 -T IndelRealigner -R $fa -I ${prefix}_rg.bam -targetIntervals indels.intervals -o ${prefix}.realign.bam
    #java -Xmx4g -jar ${workflow.projectDir}/bin/srma-0.1.15.jar I=${prefix}_rg.bam O=${prefix}.realign.bam R=${fa}
    samtools index ${prefix}.realign.bam
    """
}
