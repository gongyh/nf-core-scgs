process TNF_RPKM {
    tag "${fasta_file.baseName}"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dcvbin:933d4092ad6a07f0'

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
    python ${projectDir}/bin/dcvbin/myvae/mainfiles/calc_tnf_and_rpkm_2.py \
        -od tnf_and_rpkm \
        -fd "${fasta_file}" \
        -bam "${bam_file}" \
        --merge_samples
    """
}
