process TNF_RPKM {
    tag "${fasta_file.baseName}"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/pycoverm_pysam_loguru_numpy_pruned:10462dcc32c4ddbf'

    input:
    path fasta_file
    path bam_file

    output:
    path "tnf_and_rpkm/*tnf.npz",  emit: tnf
    path "tnf_and_rpkm/*rpkm.npz", emit: rpkm

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${fasta_file.baseName}"
    """
    mkdir -p tnf_and_rpkm
    grep ">" "${fasta_file}" | sed 's/>//g' > seq_ids.txt
    samtools view -b -h "${bam_file}" -N seq_ids.txt > filtered_perfect.bam
    samtools index filtered_perfect.bam
    python ${projectDir}/bin/dcvbin/myvae/mainfiles/calc_tnf_and_rpkm_2.py \
        -od tnf_and_rpkm \
        -fd "${fasta_file}" \
        -bam "${bam_file}"
    """
}
