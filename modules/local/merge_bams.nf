nextflow.enable.types = true

process MERGE_BAMS {
    tag "merge_bams"
    label 'process_medium'

    conda "bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    bam_files: List<Path>

    output:
    record(merged_bam: file('merged.bam'), merged_bai: file('merged.bam.bai'))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    """
    samtools merge -@ ${task.cpus} -r merged.bam ${bam_files.join(' ')}
    samtools index merged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
