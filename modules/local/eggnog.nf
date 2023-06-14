process EGGNOG {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(faa)
    path db

    output:
    path("${prefix}.emapper.annotations")

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    set +u
    source activate base
    emapper.py -i $faa -o $prefix --data_dir $db --dmnd_db $db/eggnog_proteins.dmnd -m diamond --cpu ${task.cpus}
    """
}
