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
    grep '^>' "${fasta}" | sed 's/^>//' | awk '{print \$1}' > ${prefix}_keep_names.txt
    awk '{print \$1 "\t0\t1"}' ${prefix}_keep_names.txt > ${prefix}_contig_names.bed
    samtools view -H "${bam}" > ${prefix}_full_header.sam
    grep -v '^@SQ' ${prefix}_full_header.sam > ${prefix}_header_no_sq.sam
    samtools faidx "${fasta}"
    awk '{print "@SQ\\tSN:"\$1"\\tLN:"\$2}' "${fasta}.fai" > ${prefix}_new_sq.sam
    cat ${prefix}_header_no_sq.sam ${prefix}_new_sq.sam > ${prefix}_new_header.sam
    if ! grep -q '^@SQ' ${prefix}_new_header.sam; then
        echo "ERROR: new_header.sam missing @SQ" >&2
        exit 1
    fi
    samtools index "${bam}"
    samtools view -b -L ${prefix}_contig_names.bed -o ${prefix}_body.bam "${bam}"
    samtools reheader ${prefix}_new_header.sam ${prefix}_body.bam > ${prefix}.bam
    samtools index ${prefix}.bam
    rm ${prefix}_keep_names.txt ${prefix}_contig_names.bed ${prefix}_full_header.sam ${prefix}_header_no_sq.sam ${prefix}_new_sq.sam ${prefix}_new_header.sam ${prefix}_body.bam "${fasta}.fai"




    cat <<EOF > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed 's/samtools //')
    EOF
    """
}
