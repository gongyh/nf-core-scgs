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
    record(results: file("quast_*"), tsv: file("quast_*/*.tsv"), versions: file("versions.yml"))

    script:
    def euk_cmd = euk ? (fungus ? "--fungus" : "-e") : ""
    def ref = fasta.exists() ? "-r $fasta" : ""
    def gene = gff.exists() ? "--features gene:$gff" : ""
    def outdir = quast_outdir.replaceAll(/[\\/:*?"<>|]/, '_').replaceAll(/[\s_]+/, '_').trim()
    """
    bams=($bam)
    bams_param=\$(echo \${bams[*]} | sed 's/ /,/g')
    labels=\$(echo \${bams[*]} | sed 's/.markdup.bam//g' | sed 's/ /,/g')
    quast.py -o $outdir $ref $gene -m 200 -t ${task.cpus} $euk_cmd --rna-finding --bam \$bams_param -l \$labels --no-sv --no-read-stats $contigs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
