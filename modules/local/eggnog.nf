process EGGNOG {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::eggnog-mapper=2.1.11=pyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.8--pyhdfd78af_0' :
        'biocontainers/eggnog-mapper:2.1.11--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(faa)
    path db

    output:
    tuple val(meta), path("${prefix}.emapper.annotations"), emit: annotations
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    emapper.py -i $faa -o $prefix --data_dir $db --dmnd_db $db/eggnog_proteins.dmnd -m diamond --cpu ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$(echo \$(emapper.py --version | grep emapper 2>&1 ) | cut -d'/' -f1 | sed 's/^.*emapper-//; s/Using.*\$//')
    END_VERSIONS
    """
}
