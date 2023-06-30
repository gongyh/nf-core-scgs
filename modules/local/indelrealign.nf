process INDELREALIGN {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::gatk=3.8=hdfd78af_11 bioconda::samtools=1.17 bioconda::picard=2.19.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1078cca9be35f4e979036df5cc474bba6aa9e104:a682d8c0dfb0582e13b281ab918a1de0b7a7778d-0' :
        'scgs/mulled-v2-1078cca9be35f4e979036df5cc474bba6aa9e104:a682d8c0dfb0582e13b281ab918a1de0b7a7778d-0' }"

    input:
    tuple val(meta), path(bam)
    path fa

    output:
    tuple val(meta), path("*.realign.bam"),         emit: bam
    tuple val(meta), path("*.realign.bam.bai"),     emit: bai
    path "versions.yml",                            emit: versions

    when:
    params.snv && !params.nanopore

    script:
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1 | grep SNAPSHOT | cut -d'-' -f1))
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gatk3: \$(echo \$(gatk3 -version 2>&1)
    END_VERSIONS
    """
}
