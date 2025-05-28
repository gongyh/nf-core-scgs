process METACOMPASS {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::metacompass=1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metacompass:1.12' :
        'scgs/metacompass:1.12' }"

    input:
    tuple val(meta), path(reads)
    path(ref_fna)

    output:
    tuple val(meta), path("${prefix}.metacompass.ctg.fa")                  , emit: contig
    tuple val(meta), path("${prefix}.metacompass_out")                     , emit: assembly
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def rcl = meta.single_end ? "-U ${reads[0]}" : "-P ${reads[0]},${reads[1]}"
    """
    go_metacompass.py -r ${ref_fna} ${rcl} -m 1 -g 100 -e ${prefix} -t ${task.cpus} -o ${prefix}.metacompass_out
    ln -s ${prefix}.metacompass_out/metacompass.final.ctg.fa ${prefix}.metacompass.ctg.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaCompass: 1.12
    END_VERSIONS
    """
}
