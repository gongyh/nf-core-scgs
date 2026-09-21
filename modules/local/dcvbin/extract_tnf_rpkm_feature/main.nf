nextflow.enable.types = true

process TNF_RPKM {
    tag "$meta.id"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dcvbin:933d4092ad6a07f0'

    input:
    tuple(meta: Map, fasta_file: Path, bam_file: Path)

    output:
    record(meta: meta, tnf: file('tnf_and_rpkm/*tnf.npz'), rpkm: file('tnf_and_rpkm/*rpkm.npz'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    def args    = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """

    mkdir -p tnf_and_rpkm
    python ${projectDir}/bin/dcvbin/myvae/mainfiles/calc_tnf_and_rpkm_2.py \
        -od tnf_and_rpkm \
        -fd "${fasta_file}" \
        -bam "${bam_file}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dcvbin: 732ee4257f7da200994f4c105e9dacbc74242883
    END_VERSIONS
    """
}
