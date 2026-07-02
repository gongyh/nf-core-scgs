process FEATURE_FUSION {
    tag "$meta.id"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dcvbin:933d4092ad6a07f0'

    input:
    tuple val(meta), path(fpf_file)
    tuple val(meta), path(tnf_file)
    tuple val(meta), path(rpkm_file)

    output:
    tuple val(meta), path("${prefix}_vae_features.npy"), emit: features
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/dcvbin/myvae/mainfiles/vaeTest_2.py \
        -dd "${fpf_file}" \
        -td "${tnf_file}" \
        -rd "${rpkm_file}" \
        -vd "${prefix}_vae_features.npy"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dcvbin: \$(python -c "import dcvbin; print(dcvbin.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
