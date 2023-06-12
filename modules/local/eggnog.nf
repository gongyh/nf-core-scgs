process EGGNOG {
    tag "$prefix"
    publishDir "${params.outdir}/eggnog", mode: 'copy'

    input:
    path faa
    path db

    output:
    path("${prefix}.emapper.annotations")

    when:
    eggnog_db

    script:
    prefix = faa.toString() - ~/(\.proteins\.fa)?(\.faa)?$/
    """
    set +u
    source activate base
    emapper.py -i $faa -o $prefix --data_dir $db --dmnd_db $db/eggnog_proteins.dmnd -m diamond --cpu ${task.cpus}
    """
}
