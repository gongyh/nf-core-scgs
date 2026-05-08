process EXTRACT_BINS {
    tag "extract_bins"
    label 'process_low'
    publishDir "${params.outdir}/extracted_bins", mode: 'copy'

    conda "bioconda::seqtk=1.3"
    container "community.wave.seqera.io/library/seqtk:r93--b54ec2a2e8839010"

    input:
    path clusters
    path assembly

    output:
    path "bins", emit: bins
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p bins
    if [ ! -f clusters.tsv ]; then
        echo "ERROR: clusters.tsv not found!" >&2
        exit 1
    fi
    awk 'NR>1 && \$2 != "unbinned" {print \$2}' clusters.tsv | sort -u > bin_names.txt
    if [ ! -s bin_names.txt ]; then
        echo "No bins to extract (all unbinned). Exiting gracefully." >&2
        exit 0
    fi
    while read bin; do
        echo "Processing bin: \$bin"
        awk -v b="\$bin" '\$2 == b {print \$1}' clusters.tsv > \${bin}_list.txt
        seqtk subseq ${assembly} \${bin}_list.txt > bins/\${bin}.fa
        rm \${bin}_list.txt
    done < bin_names.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | grep -oP 'Version \\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}
