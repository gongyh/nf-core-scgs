// modules/local/cooccurrence_binning.nf
// 功能：基于覆盖度矩阵进行共现分箱（Fisher检验 + t-SNE + DBSCAN）
// 输入：coverage_matrix.tsv, super_contigs.fasta
// 输出：clusters.tsv, bins/ 目录

process COOCCURRENCE_BINNING {
    tag "binning"
    label 'process_high'          // 需要较多内存和 CPU
    publishDir "${params.outdir}/binning", mode: 'copy'

    input:
    path coverage_matrix
    path super_contigs

    output:
    path "clusters.tsv", emit: clusters
    path "bins/", emit: bins
    path "versions.yml", emit: versions

    script:
    def threshold = params.coverage_thresh ?: 2048
    def eps       = params.dbscan_eps      ?: 2.6
    def minpts    = params.dbscan_minpts   ?: 5
    def seed      = params.random_seed     ?: 2015
    """
    cooccurrence_binning.py \\
        --coverage ${coverage_matrix} \\
        --contigs ${super_contigs} \\
        --outdir . \\
        --threshold ${threshold} \\
        --eps ${eps} \\
        --minpts ${minpts} \\
        --random_seed ${seed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        script: "cooccurrence_binning.py"
        numpy: \$(python -c "import numpy; print(numpy.__version__)" 2>/dev/null || echo 'N/A')
        scipy: \$(python -c "import scipy; print(scipy.__version__)" 2>/dev/null || echo 'N/A')
        sklearn: \$(python -c "import sklearn; print(sklearn.__version__)" 2>/dev/null || echo 'N/A')
        biopython: \$(python -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo 'N/A')
    END_VERSIONS
    """
}
