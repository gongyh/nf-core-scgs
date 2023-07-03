process QUAST_REF {
    label 'QUAST'

    conda "bioconda::quast=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
    path(fasta)
    path(gff)
    path(contig)
    path(bam)
    path(bai)
    val(euk)
    val(fungus)

    output:
    path "quast"       , emit: results
    path 'quast/*.tsv' , emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def euk_cmd = euk ? ( params.fungus ? "--fungus" : "-e") : ""
    def ref = fasta.exists() ? "-r $fasta" : ""
    def gene = gff.exists() ? "--features gene:$gff" : ""
    """
    contigs=\$(ls *.ctgs.fasta | paste -sd " " -)
    labels=\$(ls *.ctgs.fasta | paste -sd "," - | sed 's/.ctgs.fasta//g')
    bams=\$(ls *.markdup.bam | paste -sd "," -)
    quast.py -o quast $ref $gene -m 200 -t ${task.cpus} $euk_cmd --rna-finding --bam \$bams -l \$labels --no-sv --no-read-stats \$contigs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
