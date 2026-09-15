process FILTER_CONTIGS {
    tag "filter_contigs"
    label 'process_low'

    conda "bioconda::seqkit=2.3.1"
    container "quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0"

    input:
    path fasta
    val min_len

    output:
    path "filtered.fasta", emit: filtered
    path "versions.yml", emit: versions
    script:
    """
    seqkit seq -m ${min_len} ${fasta} > filtered.fasta
    echo "seqtk: \$(seqtk 2>&1 | head -1)" > versions.yml
    """
}
