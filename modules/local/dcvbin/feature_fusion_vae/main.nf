process FEATURE_FUSION {
    tag "${fpf_file.baseName}"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dcvbin:ea1d53670b689bf9'

    input:
    path fpf_file
    path tnf_file
    path rpkm_file

    output:
    path "${prefix}_vae_features.npy", emit: features

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${fpf_file.baseName.replaceAll('_fpf', '')}"
    """
    python ${projectDir}/bin/dcvbin/myvae/mainfiles/vaeTest_2.py \
        -dd "${fpf_file}" \
        -td "${tnf_file}" \
        -rd "${rpkm_file}" \
        -vd "${prefix}_vae_features.npy"
    """
}
