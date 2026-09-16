process MMSEQS2SEMIBIN {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("${prefix}_semibin_tax.tsv"), emit: tax
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk -F'\t' 'NR>1 {print \$1"\t"\$3}' ${tsv} > ${prefix}_semibin_tax.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -1)
    END_VERSIONS
    """
}
