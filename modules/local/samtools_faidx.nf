nextflow.enable.types = true

process SAMTOOLS_FAIDX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.23.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"
    input:
    tuple(meta: Map, fasta: Path)

    output:
    record(meta: meta, fai: file('*.fai'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    samtools faidx $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
