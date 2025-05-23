process SPADES {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::spades=3.15.5 bioconda::perl-bioperl=1.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-23d30bd0f79edd4339b884a2320935a5a236f7eb:824e273bd5969e5d2f8d617c66ab71e506b4ea71-0' :
        'scgs/mulled-v2-23d30bd0f79edd4339b884a2320935a5a236f7eb:824e273bd5969e5d2f8d617c66ab71e506b4ea71-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.contigs.fasta")                       , emit: contig
    tuple val(meta), path("${prefix}.correct.paths")                       , emit: contig_path
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
    mv ${prefix}.spades_out/assembly_graph_after_simplification.gfa ${prefix}.spades_out/${prefix}.contigs.gfa
    ln -s ${prefix}.spades_out/contigs.paths ${prefix}.contigs.paths
    ln -s ${prefix}.spades_out/contigs.fasta ${prefix}.contigs.fasta
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/_length.*\$//g' > ${prefix}.ctgs.fasta
    correctPaths.py ${prefix}.contigs.fasta ${prefix}.contigs.paths

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//; s/Using.*\$//')
    END_VERSIONS
    """
}
