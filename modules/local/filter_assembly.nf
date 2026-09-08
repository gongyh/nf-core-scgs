process FILTER_ASSEMBLY {
    tag "filter_assembly"
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
    
    cat <<-END_VERSIONS > versions.yml
    "NFCORE_MINIMETA:MINIMETA:FILTER_ASSEMBLY":
        seqkit: \$(seqkit version 2>&1 | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
