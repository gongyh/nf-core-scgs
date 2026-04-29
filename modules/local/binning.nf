process COOCCURRENCE_BINNING {
    tag "cooccurrence"
    label 'process_low'

    conda "conda-forge::python=3.9 conda-forge::pandas conda-forge::scipy conda-forge::scikit-learn bioconda::samtools"
    container "community.wave.seqera.io/library/samtools_pandas_scikit-learn_scipy:01d109078e7b563b"

    input:
    path coverage_tsv

    output:
    path "clusters.tsv", emit: clusters
    path "versions.yml", emit: versions

    script:
    def script_path = "${projectDir}/bin/cooccurrence_binning.py"
    def eps = task.ext.eps ?: 0.05
    def min_samples = task.ext.min_samples ?: 2
    """
    python ${script_path} ${coverage_tsv} clusters.tsv --eps ${eps} --min_samples ${min_samples}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
        pandas: \$(python -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo "N/A")
        scipy: \$(python -c "import scipy; print(scipy.__version__)" 2>/dev/null || echo "N/A")
        sklearn: \$(python -c "import sklearn; print(sklearn.__version__)" 2>/dev/null || echo "N/A")
        samtools: \$(samtools --version 2>&1 | head -1)
    END_VERSIONS
    """
}
