nextflow.enable.types = true

process CONTIG_KMER {
    tag "$meta.id"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dcvbin:933d4092ad6a07f0'

    input:
    tuple(meta: Map, fasta_file: Path)

    output:
    record(meta: meta, kmer: file('*4mer.csv'), seqid: file('*seqid.csv'), versions: file('versions.yml'))

    script:
    def args    = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/dcvbin/scripts/calculate_kmer_multi_thread_2.py \
        "${fasta_file}" \
        "4mer.csv" \
        "seqid.csv" \
        4

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dcvbin: \$(python -c "import dcvbin; print(dcvbin.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
