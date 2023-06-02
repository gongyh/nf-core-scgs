process SAVE_REFERENCE {
    publishDir path: "${params.outdir}/reference_genome", mode: 'copy'

    input:
    path fasta
    path gff

    output:
    path("genome.fa")
    path("genome.gff")
    path("*.bed")
    path("genome.bed"),                     emit: bed

    when:
    params.fasta && params.gff

    script:
    """
    ln -s ${fasta} genome.fa
    ln -s ${gff} genome.gff
    fa2bed.py genome.fa
    cat genome.gff | grep \$'\tgene\t' | bedtools sort | cut -f1,4,5,7 > genes.bed
    """
}
