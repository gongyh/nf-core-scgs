process CONTIG_KMER {
    tag "${fasta_file.baseName}"

    conda "${moduleDir}/dcvbin.yaml"
    container 'community.wave.seqera.io/library/dnaberts:7a7299083f265248'

    input:
    path fasta_file

    output:
    path "4mer.csv",  emit: kmer
    path "seqid.csv"

    script:
    def args    = task.ext.args ?: ''
    """
    python ${projectDir}/bin/dcvbin/scripts/calculate_kmer_multi_thread_2.py \
        "${fasta_file}" \
        "4mer.csv" \
        "seqid.csv" \
        4
    """
}
