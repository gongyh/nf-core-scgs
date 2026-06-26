process CONTIG_EMBEDDING {
    tag "${task.ext.prefix ?: ctgs_2k.baseName}"
    label 'process_gpu'

    conda "${moduleDir}/dnaberts.yaml"

    container 'community.wave.seqera.io/library/dcvbin:933d4092ad6a07f0'

    input:
    tuple val(meta), path(ctgs_2k)
    path model_dir

    output:
    path "${task.ext.prefix ?: ctgs_2k.baseName}_fpf.npy", emit: fpf

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${ctgs_2k.baseName}"
    """
    export QT_XCB_GL_INTEGRATION="none"
    python ${projectDir}/bin/dcvbin/scripts/featureExtract_gpu_2.py \
        -md "${model_dir}" \
        -fd "${ctgs_2k}" \
        -sd "${prefix}_over2kseq.txt" \
        -dd "${prefix}_fpf.npy"
    """
}
