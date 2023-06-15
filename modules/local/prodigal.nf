process PRODIGAL {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::prodigal=2.6.3 conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2e442ba7b07bfa102b9cf8fac6221263cd746ab8:57f05cfa73f769d6ed6d54144cb3aa2a6a6b17e0-0' :
        'biocontainers/mulled-v2-2e442ba7b07bfa102b9cf8fac6221263cd746ab8:57f05cfa73f769d6ed6d54144cb3aa2a6a6b17e0-0' }"

    input:
    tuple val(meta), path(contigs)

    output:
    path("$prefix")
    path "versions.yml",  emit: versions

    when:
    !euk

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    prodigal -i $contigs -o ${prefix}/${prefix}.gbk -a ${prefix}/${prefix}.proteins.faa -p meta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
    END_VERSIONS
    """
}
