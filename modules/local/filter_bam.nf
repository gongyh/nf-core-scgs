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
    samtools faidx "${fasta}"
    samtools view -H "${bam}" | grep -E "^@HD|^@PG|^@CO" > new_header.sam
    awk '{print "@SQ\\tSN:"\$1"\\tLN:"\$2}' "${fasta}.fai" >> new_header.sam
    grep '^>' "${fasta}" | sed 's/^>//' | awk '{print \$1}' > keep_names.txt
    samtools index "${bam}"
    samtools view -b -h -N keep_names.txt -o body.bam "${bam}"
    samtools reheader new_header.sam body.bam > ${prefix}.bam
    samtools index ${prefix}.bam
    rm new_header.sam body.bam keep_names.txt "${fasta}.fai"

    cat <<EOF > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    EOF
    """
}
