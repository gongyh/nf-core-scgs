process MAKE_DUMMY_TAXONOMY {
    tag "dummy_taxonomy"
    label 'process_single'

    input:
    path fasta

    output:
    path "dummy_taxonomy.tsv", emit: taxonomy

    script:
    """
    echo -e "contigs\\tpredictions" > dummy_taxonomy.tsv
    grep "^>" ${fasta} | sed 's/^>//' | awk '{print \$1"\t1"}' >> dummy_taxonomy.tsv
    """
}
