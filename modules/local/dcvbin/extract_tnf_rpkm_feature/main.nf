process TNF_RPKM {
    tag "$meta.id"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dcvbin:933d4092ad6a07f0'

    input:
    tuple val(meta), path(fasta_file)
    path bam_file

    output:
    tuple val(meta), path("tnf_and_rpkm/*tnf.npz"), emit: tnf
    tuple val(meta), path("tnf_and_rpkm/*rpkm.npz"), emit: rpkm
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

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
        dcvbin: \$(python -c "import dcvbin; print(dcvbin.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
