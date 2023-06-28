process BLASTN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::blast=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blast:2.13.0--hf3cf87c_0' :
        'biocontainers/blast:2.13.0--hf3cf87c_0' }"

    input:
    tuple val(meta), path(contigs)
    path db
    val(evalue)

    output:
    tuple val(meta), path("${contigs}" ),       emit: contigs
    tuple val(meta), path("${prefix}_nt.out"),  emit: nt
    path "versions.yml",                        emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export BLASTDB=$db
    blastn -query $contigs -db $db/nt -outfmt '6 qseqid staxids bitscore std' \
        -max_target_seqs 1 -max_hsps 1 -evalue ${evalue} \
        -num_threads ${task.cpus} -out ${prefix}_nt.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blastn: \$(echo \$(blastn -version 2>&1) | sed 's/^.*blastn //; s/Using.*\$//')
    END_VERSIONS
    """
}
