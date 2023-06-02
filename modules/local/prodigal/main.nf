process PRODIGAL {
    tag "$prefix"
    label 'process_single'
    publishDir "${params.outdir}/prodigal", mode: 'copy'

    conda "bioconda::prodigal=2.6.3 conda-forge::pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-2e442ba7b07bfa102b9cf8fac6221263cd746ab8:57f05cfa73f769d6ed6d54144cb3aa2a6a6b17e0-0' :
        'biocontainers/mulled-v2-2e442ba7b07bfa102b9cf8fac6221263cd746ab8:57f05cfa73f769d6ed6d54144cb3aa2a6a6b17e0-0' }"

    input:
    path contigs

    output:
    path("$prefix")

    when:
    !euk

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    """
    mkdir -p ${prefix}
    prodigal -i $contigs -o ${prefix}/${prefix}.gbk -a ${prefix}/${prefix}.proteins.faa -p meta
    """
}
