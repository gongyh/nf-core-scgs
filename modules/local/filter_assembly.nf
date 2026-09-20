nextflow.enable.types = true

process FILTER_ASSEMBLY {
    tag "filter_assembly"
    label 'process_low'

    conda "bioconda::seqkit=2.3.1"
    container "quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0"

    input:
    fasta: Path
    min_len: Integer

    output:
    record(filtered: file('filtered.fasta'), filtered_ids: file('filtered_ids.txt'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    seqkit seq -m ${min_len} ${fasta} > filtered.fasta
    seqkit seq -m ${min_len} ${fasta} -n -i > filtered_ids.txt
    cat <<-END_VERSIONS > versions.yml
    "NFCORE_MINIMETA:MINIMETA:FILTER_ASSEMBLY":
        seqkit: \$(seqkit version | sed 's/seqkit //')
    END_VERSIONS
    """
}
