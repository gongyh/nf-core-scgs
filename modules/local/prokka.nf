process PROKKA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::prokka=1.14.6 bioconda::bedops=2.4.38"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1e40df84b5b2d0a934c357a759500c269d2eb793:81460e1910925aa1427c823417f44d2739507564-0' :
        'scgs/mulled-v2-1e40df84b5b2d0a934c357a759500c269d2eb793:81460e1910925aa1427c823417f44d2739507564-0' }"

    input:
    tuple val(meta), path(contigs)
    path proteins

    output:
    tuple val(meta), path("$prefix")                , emit: prokka_for_split
    tuple val(meta), path("${prefix}/${prefix}.faa"), emit: faa
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    """
    cat $contigs | sed 's/_length.*\$//g' > ${prefix}_node.fa
    prokka --outdir $prefix --prefix $prefix --strain $prefix --addgenes --addmrna --cpus ${task.cpus} $proteins_opt ${prefix}_node.fa
    sed '/^##FASTA/Q' ${prefix}/${prefix}.gff > ${prefix}/${prefix}_noseq.gff
    gff2bed < ${prefix}/${prefix}_noseq.gff | cut -f1,4 | grep _gene | sed 's/_gene//g' > ${prefix}/${prefix}_ctg_genes.tsv
    prokka_postprocess.py ${prefix}/${prefix}_ctg_genes.tsv ${prefix}/${prefix}.tsv > ${prefix}/${prefix}_all.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka -v 2>&1) | sed 's/^.*prokka //; s/Using.*\$//')
    END_VERSIONS
    """
}
