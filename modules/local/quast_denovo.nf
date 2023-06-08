process QUAST_DENOVO {
    label 'process_medium'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path("*.ctgs.fasta")
    val(euk)

    output:
    path("quast/report.tsv"),  emit: report
    path("quast")

    when:
    denovo == true

    script:
    def euk_cmd = euk ? ( params.fungus ? "--fungus" : "-e") : ""
    """
    contigs=\$(ls *.ctgs.fasta | paste -sd " " -)
    labels=\$(ls *.ctgs.fasta | paste -sd "," - | sed 's/.ctgs.fasta//g')
    quast.py -o quast -m 200 -t ${task.cpus} $euk_cmd --rna-finding -l \$labels --no-sv --no-read-stats \$contigs
    """
}
