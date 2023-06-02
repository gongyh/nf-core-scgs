process BLASTN {
    tag "$prefix"
    label 'process_medium'
    publishDir "${params.outdir}/blob", mode: 'copy'

    conda "bioconda::blast=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.13.0--hf3cf87c_0' :
        'biocontainers/blast:2.13.0--hf3cf87c_0' }"

    input:
    path contigs
    path nt

    output:
    path("${prefix}_nt.out"),                emit: nt
    path("${contigs}"),                      emit: contigs

    when:
    params.nt_db

    script:
    prefix = contigs.toString() - ~/(\.ctg200\.fasta)?(\.ctg200)?(\.fasta)?(\.fa)?$/
    """
    export BLASTDB=$nt
    blastn -query $contigs -db $nt/nt -outfmt '6 qseqid staxids bitscore std' \
    -max_target_seqs 1 -max_hsps 1 -evalue ${params.evalue} \
    -num_threads ${task.cpus} -out ${prefix}_nt.out
    """
}
