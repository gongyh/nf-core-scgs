process GRAPHBIN {
    label 'process_medium'

    conda "bioconda::graphbin=1.7.1--pyh7cba7a3_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphbin:1.7.1--pyh7cba7a3_0' :
        'biocontainers/graphbin:1.7.1--pyh7cba7a3_0' }"

    input:
    path("contigs/*")
    path("paths/*")
    path("gfa/*")
    path("csv/*")

    output:
    path("binning/*")      , emit: out_put
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir binning
    samples=(`ls contigs/*.contigs.fasta | sed 's~contigs/~~g;s~.contigs.fasta~~g'`)
    for sample in \${samples[*]}; do
        mkdir binning/\$sample
        graphbin --assembler spades --graph gfa/\$sample.contigs.gfa --contigs contigs/\$sample.contigs.fasta \\
        --paths paths/\$sample.correct.paths --binned csv/\$sample.bin.csv --output binning/\$sample
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphbin: \$(echo \$(graphbin -v 2>&1) | sed 's/^.*graphbin, version //; s/Using.*\$//')
    END_VERSIONS
    """
}
