process PROKKA {
    tag "$prefix"
    label 'process_low'
    publishDir "${params.outdir}/prokka", mode: 'copy'

    conda "bioconda::prokka=1.14.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prokka%3A1.14.6--pl5321hdfd78af_4' :
        'biocontainers/prokka:1.14.6--pl5321hdfd78af_4' }"

    input:
    path contigs

    output:
    path("$prefix"),                        emit: prokka_for_split
    path("$prefix/${prefix}.faa"),          emit: faa

    when:
    !euk

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    """
    if [ \$(id -u) -eq 0 ]; then
    wget -c -q -t 1 -T 60 ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz -O linux64.tbl2asn.gz && gunzip linux64.tbl2asn.gz && chmod +x linux64.tbl2asn && mv linux64.tbl2asn /opt/conda/envs/nf-core-gongyh-scgs/bin/tbl2asn
    fi
    cat $contigs | sed 's/_length.*\$//g' > ${prefix}_node.fa
    prokka --outdir $prefix --prefix $prefix --addgenes --cpus ${task.cpus} ${prefix}_node.fa || echo "Ignore minor errors of prokka!"
    sed '/^##FASTA/Q' ${prefix}/${prefix}.gff > ${prefix}/${prefix}_noseq.gff
    gff2bed < ${prefix}/${prefix}_noseq.gff | cut -f1,4 | grep -v gene > ${prefix}/${prefix}_ctg_genes.tsv
    prokka_postprocess.py ${prefix}/${prefix}_ctg_genes.tsv ${prefix}/${prefix}.tsv > ${prefix}/${prefix}_all.tsv
    """
}
