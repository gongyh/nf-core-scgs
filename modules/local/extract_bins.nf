process EXTRACT_BINS {
    tag "extract_bins"
    label 'process_low'

    conda "bioconda::seqtk=1.3"
    container "community.wave.seqera.io/library/seqtk:r93--b54ec2a2e8839010"

    input:
    path clusters
    path assembly

    output:
    path "bins", emit: bins
    path "versions.yml", emit: versions
    path "extract_bins_mqc.tsv", emit: mqc_tsv
    path "scaffolds2bin.tsv", emit: scaffolds2bin
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

    > scaffolds2bin.tsv
    for bin_fa in bins/*.fa; do
        if [ -f "\$bin_fa" ]; then
            bin_name=\$(basename "\$bin_fa" .fa)
            grep "^>" "\$bin_fa" | sed 's/^>//' | awk -v bin="\$bin_name" '{print \$1"\t"bin}'
        fi
    done >> scaffolds2bin.tsv
    if [ -d "bins" ]; then
        N_BINS=\$(ls bins/*.fa 2>/dev/null | wc -l)
        TOTAL_SIZE=\$(ls -l bins/*.fa 2>/dev/null | awk '{sum+=\$5} END {print sum}')
        [ -z "\$TOTAL_SIZE" ] && TOTAL_SIZE=0
    else
        N_BINS=0; TOTAL_SIZE=0
    fi

    printf "Metric\tValue\n" > extract_bins_mqc.tsv
    printf "Number of bins extracted\t\${N_BINS}\n" >> extract_bins_mqc.tsv
    printf "Total bin size (bp)\t\${TOTAL_SIZE}\n" >> extract_bins_mqc.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | grep -oP 'Version \\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """
}
