process CLASSIFY_TAXA {
    tag "classify_taxa"
    label 'process_single'

    conda "conda-forge::coreutils=9.4"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/alpine:3.19' :
        'alpine:3.19' }"
    input:
    path fasta

    output:
    path "taxonomy.tsv", emit: taxonomy
    path "versions.yml", emit: versions

    script:
    """
    echo -e "contigs\\tpredictions\\tscores" > taxonomy.tsv
    grep "^>" ${fasta} | sed 's/^>//' | awk '{print \$1"\t1\t1.0"}' >> taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        classify_taxa: \$(echo "dummy")
    END_VERSIONS
    """
}
