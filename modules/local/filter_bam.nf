process FILTER_BAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::samtools=1.21"
    container 'community.wave.seqera.io/library/dcvbin:ea1d53670b689bf9'

    input:
    tuple val(meta), path(fasta)
    path(bam)
    path(bai)

    output:
    tuple val(meta), path("${meta.id}_filtered.bam"), emit: filtered_bam
    path "${meta.id}_filtered.bam.bai", emit: filtered_bai
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}_filtered"
    """
    # 1. 动态为过滤后的 fasta 文件创建 FAI 索引
    samtools faidx "${fasta}"

    # 2. 【核心保底】：利用 samtools mawk/awk 或者内置命令，直接用新 fasta 的索引动态生成一个绝对纯净的新 Header 文本
    # 它会动态生成只有 80 条序列的 @SQ 标签
    samtools view -H "${bam}" | grep -E "^@HD|^@PG|^@CO" > new_header.sam
    awk '{print "@SQ\\tSN:"\$1"\\tLN:"\$2}' "${fasta}.fai" >> new_header.sam

    # 3. 提取原 BAM 中属于这 80 条序列的比对 Read（临时文件）
    grep '^>' "${fasta}" | sed 's/^>//' | awk '{print \$1}' > keep_names.txt
    samtools index "${bam}"
    samtools view -b -h -N keep_names.txt -o body.bam "${bam}"

    # 4. 【灵魂合体】：用 samtools reheader 将我们纯手工、动态创建的 80 条干净 Header 强行写进 BAM
    samtools reheader new_header.sam body.bam > ${prefix}.bam

    # 5. 为最终完美的 BAM 构建新索引
    samtools index ${prefix}.bam

    # 清理临时文件
    rm new_header.sam body.bam keep_names.txt "${fasta}.fai"

    cat <<EOF > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    EOF
    """
}
