nextflow.enable.types = true

process SEMIBIN2 {
    tag "semibin2_binning"
    label 'process_medium'
    label 'process_gpu'

    conda "semibin_env.yaml"
    container 'community.wave.seqera.io/library/python_pip_bedtools_hmmer_pruned:b5edf268fd86239c'

    input:
    assembly: Path
    merged_bam: Path
    taxonomy: Path?
    output:
    record(bins: file('bins_merged'), mqc_tsv: file('semibin2_mqc.tsv'), scaffolds2bin: file('scaffolds2bin.tsv'))
    topic:
    file('versions.yml') >> 'versions'
    script:
    def bam_args = "-b ${merged_bam}"
    def tax_args = taxonomy ? "--semi-supervised --taxonomy-annotation-table ${taxonomy}" : "--self-supervised"
    def args = task.ext.args ?: ''
    """
    # Precomputed annotations bypass classification, including its dependency check.
    cat > semibin_launcher.py <<'PY'
import sys
from SemiBin import main

if __name__ == '__main__':
    mode = sys.argv.pop(1)
    if mode == 'precomputed':
        check_install = main.check_install
        def check_precomputed(*args, **kwargs):
            kwargs['allow_missing_mmseqs2'] = True
            return check_install(*args, **kwargs)
        main.check_install = check_precomputed
    main.main2(sys.argv[1:])
PY
    python semibin_launcher.py ${taxonomy ? 'precomputed' : 'self'} single_easy_bin \\
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
            grep "^>" "\$bin_fa" | sed 's/^>//' | awk -v bin="\$bin_name" '{print \$1"\\t"bin}'
        done >> scaffolds2bin.tsv
    else
        echo "WARNING: No .fa files found in bins_merged, creating empty scaffolds2bin.tsv" >&2
        touch scaffolds2bin.tsv
    fi
    N_BINS=\$(find bins_merged -maxdepth 1 -name '*.fa' | wc -l)
    printf "Metric\\tValue\\n" > semibin2_mqc.tsv
    printf "Number of bins recovered\\t\${N_BINS}\\n" >> semibin2_mqc.tsv

    printf '${task.process}:\\n  semibin2: %s\\n' "\$(SemiBin2 --version 2>&1 | head -1)" > versions.yml
    """

    stub:
    """
    mkdir -p bins_merged
    touch scaffolds2bin.tsv
    printf 'Metric\\tValue\\nNumber of bins recovered\\t0\\n' > semibin2_mqc.tsv
    printf '${task.process}:\\n  semibin2: stub\\n' > versions.yml
    """
}
