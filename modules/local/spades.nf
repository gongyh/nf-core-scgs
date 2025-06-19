process SPADES {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::spades=3.15.5 bioconda::perl-bioperl=1.7.8 conda-forge::python=3.10.14"
    container "scgs/mulled-v2-5524a20c8f39de906b127a66052c67b51c9a9ce1:c8e22953d04dee6a4da05f7a131bbd081ad78651-0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.contigs.fasta")                       , emit: contig
    tuple val(meta), path("${prefix}.contigs.paths")                       , emit: contig_path
    tuple val(meta), path("${prefix}.spades_out/${prefix}.contigs.gfa")    , emit: contig_graph
    tuple val(meta), path("${prefix}.ctg200.fasta")                        , emit: ctg200
    tuple val(meta), path("${prefix}.ctgs.fasta")                          , emit: ctg
    tuple val(meta), path("${prefix}.spades_out")                          , emit: assembly
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def mode = params.bulk ? "--cov-cutoff auto" : "--sc"
    def rcl = meta.single_end ? "-s ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    spades.py ${rcl} ${mode} ${args} -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    cp ${prefix}.spades_out/assembly_graph_after_simplification.gfa ${prefix}.spades_out/${prefix}.contigs.gfa
    cp ${prefix}.spades_out/contigs.paths ${prefix}.spades_out/contigs.paths_raw
    cp ${prefix}.spades_out/scaffolds.paths ${prefix}.spades_out/scaffolds.paths_raw
    correctPaths.py ${prefix}.spades_out/contigs.fasta ${prefix}.spades_out/contigs.paths_raw ${prefix}.spades_out/contigs.paths
    correctPaths.py ${prefix}.spades_out/scaffolds.fasta ${prefix}.spades_out/scaffolds.paths_raw ${prefix}.spades_out/scaffolds.paths
    cp ${prefix}.spades_out/contigs.paths ${prefix}.contigs.paths
    cp ${prefix}.spades_out/contigs.fasta ${prefix}.contigs.fasta
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/_length.*\$//g' > ${prefix}.ctgs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//; s/Using.*\$//')
    END_VERSIONS
    """
}
