nextflow.enable.types = true

process MARKER_NCLUSTERS {
    tag "$meta.id"

    conda "${moduleDir}/copygen.yaml"
    container 'community.wave.seqera.io/library/copygen:aca96b4a00a56131'

    input:
    tuple(meta: Map, kmer_file: Path, fasta_file: Path)

    output:
    record(meta: meta, marker_cv: file('cluster_value'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    def args    = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/dcvbin/marker_gene/src/marker_gene_utils.py \
        -kf "${kmer_file}" \
        -cf "${fasta_file}" \
        -cvf "cluster_value"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dcvbin: 732ee4257f7da200994f4c105e9dacbc74242883
    END_VERSIONS
    """
}
