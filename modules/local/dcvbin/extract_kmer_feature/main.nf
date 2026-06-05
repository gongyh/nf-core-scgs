process CONTIG_KMER {
    tag "${meta.id}"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dcvbin:ea1d53670b689bf9'

    input:
    path fasta_file

    output:
    path "4mer.csv",  emit: kmer
    path "seqid.csv"

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    python ${projectDir}/bin/dcvbin/scripts/calculate_kmer_multi_thread_2.py \
        "${fasta_file}" \
        "4mer.csv" \
        "seqid.csv" \
        4
    """
}
