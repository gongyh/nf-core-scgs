nextflow.enable.types = true

process FEATURE_FUSION {
    tag "$meta.id"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dcvbin:933d4092ad6a07f0'


    input:
    tuple(meta: Map, fpf_file: Path, tnf_file: Path, rpkm_file: Path)

    output:
    record(meta: meta, features: file("${prefix}_vae_features.npy"))
    topic:
    file('versions.yml') >> 'versions'

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
        dcvbin: 732ee4257f7da200994f4c105e9dacbc74242883
    END_VERSIONS
    """
}
