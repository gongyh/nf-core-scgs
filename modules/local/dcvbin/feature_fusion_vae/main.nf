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
    file('versions.yml') >> 'local_versions'

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
