process EGGNOG {
    tag "$meta.id"
    label 'process_medium'

    conda "eggnog-mapper=2.1.11--pyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.8--pyhdfd78af_0' :
        'biocontainers/eggnog-mapper:2.1.11--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(faa)
    path db

    output:
    path("${prefix}.emapper.annotations")
    path "versions.yml",  emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    emapper.py -i $faa -o $prefix --data_dir $db --dmnd_db $db/eggnog_proteins.dmnd -m diamond --cpu ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$(echo \$(pemapper.py --version | grep emapper 2>&1 ))
    END_VERSIONS
    """
}
