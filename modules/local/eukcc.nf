process EUKCC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::eukcc=2.1.0=pypyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukcc:2.1.0--pypyhdfd78af_0' :
        'biocontainers/eukcc:2.1.0--pypyhdfd78af_0' }"

    input:
    tuple val(meta), path(contig)
    path db

    output:
    path("${prefix}")
    path "versions.yml", emit: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    cat $contig | sed 's/_length.*\$//g' > ${prefix}_clean.fasta
    eukcc single --out $prefix --db $db --threads ${task.cpus} ${prefix}_clean.fasta || echo "Ignore minor errors of eukcc!"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eukcc: \$(echo \$(eukcc -v 2>&1) | sed 's/^.*EukCC version //; s/Using.*\$//')
    END_VERSIONS
    """
}
