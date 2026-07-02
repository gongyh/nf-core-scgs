process CONTIG_EMBEDDING {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/dnaberts.yaml"

    container 'community.wave.seqera.io/library/dcvbin:933d4092ad6a07f0'

    input:
    tuple val(meta), path(ctgs_2k)
    path model_dir

    output:
    tuple val(meta), path("${prefix}_fpf.npy"), emit: fpf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export QT_XCB_GL_INTEGRATION="none"
    python ${projectDir}/bin/dcvbin/scripts/featureExtract_gpu_2.py \
        -md "${model_dir}" \
        -fd "${ctgs_2k}" \
        -sd "${prefix}_over2kseq.txt" \
        -dd "${prefix}_fpf.npy"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dcvbin: \$(python -c "import dcvbin; print(dcvbin.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
