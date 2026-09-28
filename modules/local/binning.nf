nextflow.enable.types = true

process COOCCURRENCE_BINNING {
    tag "cooccurrence"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::pandas conda-forge::scipy conda-forge::scikit-learn conda-forge::matplotlib-base"
    container "community.wave.seqera.io/library/dnaberts:7a7299083f265248"

    input:
    coverage_tsv: Path
    filtered_ids: Path
    output:
    record(clusters: file('clusters.tsv'), mqc_tsv: file('cooccurrence_mqc.tsv'), N_BINS: env('N_BINS'), coverage_heatmap: file('coverage_heatmap.png'), pvalue_heatmap: file('pvalue_heatmap.png'), tsne_embedding: file('tsne_embedding.png', optional: true))
    topic:
    file('versions.yml') >> 'versions'
    script:
    def script_path = "${projectDir}/bin/cooccurrence_binning.py"
    def args = task.ext.args ?: ''
    def eps = params.cooccurrence_eps ?: 0.05
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    python ${script_path} \\
        ${coverage_tsv} ${filtered_ids} clusters.tsv \\
        --eps ${eps} ${args}

    if [ -f "clusters.tsv" ]; then
        N_BINS=\$(awk 'NR > 1 && \$2 != "unbinned" {bins[\$2] = 1} END {print length(bins)}' clusters.tsv)
    else
        N_BINS=0
    fi

    printf "Metric\\tValue\\n" > cooccurrence_mqc.tsv
    printf "Number of genome bins\\t\${N_BINS}\\n" >> cooccurrence_mqc.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
        pandas: \$(python -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "N/A")
        scipy: \$(python -c "import scipy; print(scipy.__version__)" 2>/dev/null || echo "N/A")
        sklearn: \$(python -c "import sklearn; print(sklearn.__version__)" 2>/dev/null || echo "N/A")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
