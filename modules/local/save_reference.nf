nextflow.enable.types = true

process SAVE_REFERENCE {
    label 'process_low'

    conda "conda-forge::click=8.1.3 conda-forge::biopython=1.81 bioconda::bedtools=2.31.0"
    container "scgs/mulled-v2-03f569b0930bbc8a26531ce48223cd6880134686:eeee3d8bada9c650a6eab38b1eecb7d20fe49a3a-0"

    input:
    fasta: Path
    gff: Path

    output:
    record(fa: file('genome.fa'), gff: file('genome.gff'), bed: file('genome.bed'), gc_bed: file('gc.bed'), gc_skew_bed: file('gcSkew.bed'), genes_bed: file('genes.bed'), versions: file('versions.yml'))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    """
    ln -s ${fasta} genome.fa
    ln -s ${gff} genome.gff
    fa2bed.py genome.fa
    if [ -s genome.gff ]; then
        cat genome.gff | grep \$'\\tgene\\t' | bedtools sort | cut -f1,4,5,7 > genes.bed
    else
        touch genes.bed
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
