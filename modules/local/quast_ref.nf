nextflow.enable.types = true

process QUAST_REF {
    tag "$quast_outdir"
    label 'process_medium'

    conda "bioconda::quast=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
    fasta: Path
    gff: Path
    contigs: Bag<Path>
    bam: Bag<Path>
    bai: Bag<Path>
    euk: Boolean
    fungus: Boolean
    quast_outdir: String

    output:
    record(results: file("quast_*"), tsv: file("quast_*/report.tsv"), versions: file("versions.yml"))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    def euk_cmd = euk ? (fungus ? "--fungus" : "-e") : ""
    def ref = fasta.exists() ? "-r $fasta" : ""
    def gene = gff.exists() ? "--features gene:$gff" : ""
    def outdir = quast_outdir.replaceAll(/[\\/:*?"<>|]/, '_').replaceAll(/[\s_]+/, '_').trim()
    def bam_param = bam.join(',')
    def labels = bam.collect { Path bam_file -> "$bam_file".replaceFirst(/\.markdup\.bam$/, '') }.join(',')
    def contig_files = contigs.join(' ')
    """
    quast.py -o $outdir $ref $gene -m 200 -t ${task.cpus} $euk_cmd --rna-finding --bam $bam_param -l $labels --no-sv --no-read-stats $contig_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
