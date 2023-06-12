process DIAMOND_BLASTX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::diamond=2.0.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_0' :
        'biocontainers/diamond:2.0.15--hb97b32f_0' }"

    input:
    tuple val(meta), path contigs
    tuple val(meta), path nt_out
    path uniprot
    path("uniprot.taxids")

    output:
    tuple val(meta), path("${prefix}_uniprot.taxified.out") , emit: uniprot
    tuple val(meta), path("${contigs}") ,                     emit: contigs
    tuple val(meta), path("${nt_out}") ,                      emit: nt
    val used ,                                                emit: real
    path("${prefix}_uniprot.*")

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if ( uniprot.toString().equals("/dev/null") || uniprot.toString().equals("null") ) {
    def used = false
    """
    touch ${prefix}_uniprot.out
    touch ${prefix}_uniprot.taxified.out
    """
    } else {
    def used = true
    """
    diamond blastx --query $contigs --db $uniprot -p ${task.cpus} -o ${prefix}_uniprot.out \
        --outfmt 6 --sensitive --max-target-seqs 1 --evalue ${params.evalue} -b ${params.blockSize}
    blobtools taxify -f ${prefix}_uniprot.out -m uniprot.taxids -s 0 -t 2
    """
    }
}
