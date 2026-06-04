process TNF_RPKM {
    tag "${meta.id}"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dcvbin:ea1d53670b689bf9'
    
    input:
    path fasta_file
    path bam_file
    
    output:
    path "tnf_and_rpkm/*tnf.npz",  emit: tnf
    path "tnf_and_rpkm/*rpkm.npz", emit: rpkm
    
    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/dcvbin/myvae/mainfiles/calc_tnf_and_rpkm_2.py \
        -od tnf_and_rpkm
        -fd "${fasta_file}" \
        -bam "${bam_file}"
    """
}
