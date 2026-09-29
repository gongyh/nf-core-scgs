nextflow.enable.types = true

process METABULI_TAXA {
    tag "metabuli_taxa"
    label 'process_medium'

    conda "bioconda::metabuli=1.2.0"
    container 'community.wave.seqera.io/library/metabuli:1.2.0--aade40d1e84cdec2'

    input:
    tuple(meta: Map, assembly: Path)
    db_dir: Path

    output:
    record(meta: meta, classifications: file("metabuli_out/${prefix}_job_classifications.tsv"), report: file("metabuli_out/${prefix}_job_report.tsv"))
    topic:
    file('versions.yml') >> 'versions'

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    metabuli classify \\
        ${assembly} \\
        ${db_dir} \\
        metabuli_out \\
        ${prefix}_job \\
        --threads ${task.cpus} \\
        ${args} \\
        --lineage 1
    printf '${task.process}:\\n  metabuli: %s\\n' "\$(metabuli --version 2>&1 | awk '/metabuli Version:/ {print \$3}')" > versions.yml
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p metabuli_out
    printf '#is_classified\\tname\\ttaxID\\tquery_length\\tscore\\te_value\\trank\\tlineage\\ttaxID:match_count\\n' > metabuli_out/${prefix}_job_classifications.tsv
    awk '/^>/ {sub(/^>/, ""); split(\$0, fields, /[ \\t]/); printf "0\\t%s\\t0\\t0\\t0\\t-\\t-\\t-\\t-\\n", fields[1]}' ${assembly} >> metabuli_out/${prefix}_job_classifications.tsv
    touch metabuli_out/${prefix}_job_report.tsv
    printf '${task.process}:\\n  metabuli: stub\\n' > versions.yml
    """
}
