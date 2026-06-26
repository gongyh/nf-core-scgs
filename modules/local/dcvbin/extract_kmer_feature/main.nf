process CONTIG_KMER {
    tag "${task.ext.prefix ?: fasta_file.baseName}"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dcvbin:933d4092ad6a07f0'

    input:
    path fasta_file

    output:
    path "4mer.csv",  emit: kmer
    path "seqid.csv"

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${fasta_file.baseName}"
    """
    python ${projectDir}/bin/dcvbin/scripts/calculate_kmer_multi_thread_2.py \
        "${fasta_file}" \
        "4mer.csv" \
        "seqid.csv" \
        4
    """
}
