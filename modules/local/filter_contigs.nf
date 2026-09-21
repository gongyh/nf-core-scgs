nextflow.enable.types = true

process FILTER_CONTIGS {
    tag "filter_contigs"
    label 'process_low'

    conda "bioconda::seqkit=2.3.1"
    container "quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0"

    input:
    fasta: Path
    min_len: Integer

    output:
    record(filtered: file('filtered.fasta'))
    topic:
    file('versions.yml') >> 'versions'
    script:
    """
    seqkit seq -m ${min_len} ${fasta} > filtered.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqkit version | sed 's/seqkit //')
    END_VERSIONS
    """
}
