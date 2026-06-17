process CONTIG_EMBEDDING {
    tag "${meta.id}"
    label 'process_gpu'

    conda "${moduleDir}/dnaberts.yaml"

    container 'community.wave.seqera.io/library/pandas_python_pytorch_scikit-learn_transformers:ac41f1afd77517cf'

    input:
    tuple val(meta), path(ctgs_2k)
    path model_dir

    output:
    path "${prefix}_fpf.npy", emit: fpf

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    export QT_XCB_GL_INTEGRATION="none"
    python ${projectDir}/bin/dcvbin/scripts/featureExtract_gpu_2.py \
        -md "${model_dir}" \
        -fd "${ctgs_2k}" \
        -sd "${prefix}_over2kseq.txt" \
        -dd "${prefix}_fpf.npy"
    """
}
