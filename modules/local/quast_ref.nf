process QUAST_REF {
    label "quast"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fasta
    path gff
    tuple val(meta), path("*.ctgs.fasta")
    tuple val(meta), path("*.markdup.bam")
    tuple val(meta), path("*.markdup.bam.bai")
    val(euk)

    output:
    path("quast/report.tsv"),  emit: report
    path("quast")

    when:
    denovo == false

    script:
    def euk_cmd = euk ? ( params.fungus ? "--fungus" : "-e") : ""
    def ref = fasta.exists() ? "-r $fasta" : ""
    def gene = gff.exists() ? "--features gene:$gff" : ""
    """
    contigs=\$(ls *.ctgs.fasta | paste -sd " " -)
    labels=\$(ls *.ctgs.fasta | paste -sd "," - | sed 's/.ctgs.fasta//g')
    bams=\$(ls *.markdup.bam | paste -sd "," -)
    quast.py -o quast $ref $gene -m 200 -t ${task.cpus} $euk_cmd --rna-finding --bam \$bams -l \$labels --no-sv --no-read-stats \$contigs
    """
}
