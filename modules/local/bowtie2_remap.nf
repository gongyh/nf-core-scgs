process BOWTIE2_REMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bowtie2=2.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0' :
        'biocontainers/bowtie2:2.4.4--py39hbb4e92a_0' }"

    input:
    tuple val(meta), path(contigs)

    output:
    path("${prefix}Bowtie2Index"), emit: index
    path "versions.yml",           emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}Bowtie2Index; cd ${prefix}Bowtie2Index
    ln -s ../${contigs} ${prefix}.fa
    bowtie2-build ${prefix}.fa ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
