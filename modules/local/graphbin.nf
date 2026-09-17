nextflow.enable.types = true

process GRAPHBIN {
    label 'process_medium'

    conda "bioconda::graphbin=1.7.1--pyh7cba7a3_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/graphbin:1.7.1--pyh7cba7a3_0' :
        'biocontainers/graphbin:1.7.1--pyh7cba7a3_0' }"

    input:
    contigs: Bag<Path>
    paths: Bag<Path>
    gfa: Bag<Path>
    csv: Bag<Path>

    output:
    record(out_put: file("binning/*"), versions: file("versions.yml"))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    """
    mkdir -p contigs paths gfa csv
    ln -s ${contigs} contigs/
    ln -s ${paths} paths/
    ln -s ${gfa} gfa/
    ln -s ${csv} csv/
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
