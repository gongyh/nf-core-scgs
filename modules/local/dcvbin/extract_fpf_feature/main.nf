nextflow.enable.types = true

process CONTIG_EMBEDDING {
    tag "$meta.id"
    label 'process_gpu'

    conda "${moduleDir}/dnaberts.yaml"

    container 'community.wave.seqera.io/library/dcvbin:933d4092ad6a07f0'

    input:
    tuple(meta: Map, ctgs_2k: Path, model_dir: Path)

    output:
    record(meta: meta, fpf: file("${prefix}_fpf.npy"))
    topic:
    file('versions.yml') >> 'versions'

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
        dcvbin: 732ee4257f7da200994f4c105e9dacbc74242883
    END_VERSIONS
    """
}
