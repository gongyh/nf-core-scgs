process QUAST_DENOVO {
    tag "$outdir"
    label 'process_medium'

    conda "bioconda::quast=5.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
    path(contig)
    val(euk)
    val(fungus)
    val(quast_outdir)

    output:
    path "${outdir}"             , emit: results
    path "${outdir}/*.tsv"       , emit: tsv
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def euk_cmd = euk ? ( params.fungus ? "--fungus" : "-e") : ""
    outdir = "${quast_outdir}".replaceAll(/[\\/:*?"<>|]/, '_').replaceAll(/[\s_]+/, '_').trim()
    """
    contigs=\$(ls *.fasta | paste -sd " " -)
    labels=\$(ls *.fasta | paste -sd "," - | sed 's/.fasta//g')
    quast.py -o $outdir -m 200 -t ${task.cpus} $euk_cmd --rna-finding -l \$labels --no-sv --no-read-stats \$contigs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
