nextflow.enable.types = true

process MMSEQS2SEMIBIN {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple(meta: Map, tsv: Path)

    output:
    record(meta: meta, tax: file("*_semibin_tax.tsv"), versions: file("versions.yml"))

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
