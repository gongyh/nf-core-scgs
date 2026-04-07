// modules/local/merge_coverage.nf
// 功能：将所有子样本的覆盖度文件合并为一个矩阵，并进行长度归一化
// 输入：cov_files（每个子样本的 coverage.tsv 文件列表），super_contigs.fasta
// 输出：coverage_matrix.tsv（行=contig，列=子样本，值为归一化覆盖度）

process MERGE_COVERAGE {
    tag "merge"
    label 'process_medium'
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
    path cov_files   // 每个子样本的 coverage.tsv 文件列表
    path fasta       // super_contigs.fasta

    output:
    path "coverage_matrix.tsv", emit: matrix
    path "versions.yml", emit: versions

    script:
    """
    python3 <<EOF
import sys
import pandas as pd
from Bio import SeqIO

# 读取每个 contig 的长度
contig_lengths = {}
for rec in SeqIO.parse("${fasta}", "fasta"):
    contig_lengths[rec.id] = len(rec.seq)

# 合并所有覆盖度文件
merged = None
for fname in "${cov_files}".split():
    df = pd.read_csv(fname, sep='\t', header=None, names=['contig', 'coverage'])
    df.set_index('contig', inplace=True)
    # 从文件名中提取子样本 ID（假设文件名格式为 {id}.coverage.tsv）
    import os
    sample_id = os.path.basename(fname).replace('.coverage.tsv', '')
    df.rename(columns={'coverage': sample_id}, inplace=True)
    if merged is None:
        merged = df
    else:
        merged = merged.join(df, how='outer')

merged.fillna(0, inplace=True)

# 归一化：覆盖度 / contig长度
for col in merged.columns:
    merged[col] = merged[col] / merged.index.map(contig_lengths)

# 保存矩阵
merged.to_csv('coverage_matrix.tsv', sep='\t')
EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        script: "merge_coverage.nf"
        pandas: \$(python -c "import pandas; print(pandas.__version__)" 2>/dev/null || echo 'N/A')
        biopython: \$(python -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo 'N/A')
    END_VERSIONS
    """
}
