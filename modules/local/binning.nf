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
    path "cooccurrence_mqc.tsv", emit: mqc_tsv
    script:
    def script_path = "${projectDir}/bin/cooccurrence_binning.py"
    def args = task.ext.args ?: ''
    """
    python ${script_path} ${coverage_tsv} clusters.tsv ${args}

    if [ -f "clusters.tsv" ]; then
        N_BINS=\$(tail -n +2 clusters.tsv | cut -f2 | sort -u | wc -l)
    else
        N_BINS=0
    fi

    printf "Metric\tValue\n" > cooccurrence_mqc.tsv
    printf "Number of genome bins\t\${N_BINS}\n" >> cooccurrence_mqc.tsv
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
