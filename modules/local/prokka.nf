process PROKKA {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("$prefix"),                 emit: prokka_for_split
    tuple val(meta), path("${prefix}/${prefix}.faa"), emit: faa
    path "versions.yml",                              emit: versions


    when:
    !euk

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ \$(id -u) -eq 0 ]; then
    wget -c -q -t 1 -T 60 ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz -O linux64.tbl2asn.gz && gunzip linux64.tbl2asn.gz && chmod +x linux64.tbl2asn && mv linux64.tbl2asn /opt/conda/envs/nf-core-gongyh-scgs/bin/tbl2asn
    fi
    cat $contigs | sed 's/_length.*\$//g' > ${prefix}_node.fa
    prokka --outdir $prefix --prefix $prefix --addgenes --cpus ${task.cpus} ${prefix}_node.fa || echo "Ignore minor errors of prokka!"
    sed '/^##FASTA/Q' ${prefix}/${prefix}.gff > ${prefix}/${prefix}_noseq.gff
    gff2bed < ${prefix}/${prefix}_noseq.gff | cut -f1,4 | grep -v gene > ${prefix}/${prefix}_ctg_genes.tsv
    prokka_postprocess.py ${prefix}/${prefix}_ctg_genes.tsv ${prefix}/${prefix}.tsv > ${prefix}/${prefix}_all.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka -v 2>&1) | sed 's/^.*prokka //; s/Using.*\$//')
    END_VERSIONS
    """
}
