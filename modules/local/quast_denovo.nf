nextflow.enable.types = true

process QUAST_DENOVO {
    tag "$quast_outdir"
    label 'process_medium'

    conda "bioconda::quast=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
    contig: Bag<Path>
    euk: Boolean
    fungus: Boolean
    quast_outdir: String

    output:
    record(results: file("quast_*"), tsv: file("quast_*/report.tsv"))
    topic:
    file("versions.yml") >> 'versions'

    script:
    def euk_cmd = euk ? (fungus ? "--fungus" : "-e") : ""
    def outdir = quast_outdir.replaceAll(/[\\/:*?"<>|]/, '_').replaceAll(/[\s_]+/, '_').trim()
    def label_suffix = quast_outdir == 'quast_spades' ? '_spades' : ''
    """
    contigs=\$(ls *.fasta | paste -sd " " -)
    labels=\$(ls *.fasta | sed 's/\\.fasta\$/${label_suffix}/' | paste -sd "," -)
    quast.py -o $outdir -m 200 -t ${task.cpus} $euk_cmd --rna-finding -l \$labels --no-sv --no-read-stats \$contigs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
