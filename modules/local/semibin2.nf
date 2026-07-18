process SEMIBIN2 {
    tag "coassembly_binning"
    label 'process_medium'
    conda "bioconda::semibin=2.3.0"
    container 'community.wave.seqera.io/library/semibin:2.3.0--33e3e4e2b94625ad'

    input:
    path assembly
    path merged_bam
    path taxonomy
    output:
    path "bins_merged", emit: bins
    path "semibin2_mqc.tsv", emit: mqc_tsv
    path "versions.yml", emit: versions
    path "scaffolds2bin.tsv", emit: scaffolds2bin
    script:
    def bam_args = "-b ${merged_bam}"
    def tax_args = taxonomy ? "--taxonomy ${taxonomy}" : ""
    def args = task.ext.args ?: ''
    """
    SemiBin2 single_easy_bin \\
        -i ${assembly} \\
        ${bam_args} \\
        ${tax_args} \\
        -o bins_merged \\
        --threads ${task.cpus} \\
        --compression none \\
        ${args}
    if [ -d bins_merged/output_bins ]; then
        mv bins_merged/output_bins/* bins_merged/ 2>/dev/null || true
        rmdir bins_merged/output_bins
    fi
    > scaffolds2bin.tsv
    if [ -d bins_merged ] && [ "\$(ls bins_merged/*.fa 2>/dev/null | wc -l)" -gt 0 ]; then
        for bin_fa in bins_merged/*.fa; do
            bin_name=\$(basename "\$bin_fa" .fa)
            grep "^>" "\$bin_fa" | sed 's/^>//' | awk -v bin="\$bin_name" '{print \$1"\t"bin}'
        done >> scaffolds2bin.tsv
    else
        echo "WARNING: No .fa files found in bins_merged, creating empty scaffolds2bin.tsv" >&2
        touch scaffolds2bin.tsv
    fi
    N_BINS=\$(find bins_merged -maxdepth 1 -name '*.fa' | wc -l)
    printf "Metric\tValue\\n" > semibin2_mqc.tsv
    printf "Number of bins recovered\t\${N_BINS}\\n" >> semibin2_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        semibin2: \$(SemiBin2 --version 2>&1 | head -1)
    END_VERSIONS
    """
}
