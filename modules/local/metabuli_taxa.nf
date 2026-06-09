process METABULI_TAXA {
    tag "metabuli_taxa"
    label 'process_medium'

    conda "bioconda::metabuli"
    container 'community.wave.seqera.io/library/metabuli:1.2.0--aade40d1e84cdec2'

    input:
    tuple val(meta), path(assembly)
    path  db_dir

    output:
    tuple val(meta), path("taxonomy.tsv"), emit: taxonomy
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    DB_DIR="${db_dir}"
    echo "DEBUG: Using database directory: \$DB_DIR"
    ls -la "\$DB_DIR"

    metabuli classify \
        ${assembly} \
        "\$DB_DIR" \
        metabuli_out \
        ${prefix}_job \
        --threads ${task.cpus} \
        ${args} \
        --lineage 1
    awk 'BEGIN {print "contigs\tpredictions"} NR>1 && !/^#/ {print \$2"\t"\$3}' metabuli_out/${prefix}_job_classifications.tsv > taxonomy.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metabuli: \$(metabuli --version 2>&1 | awk '/metabuli Version:/ {print \$3}')
    END_VERSIONS
    """
}
