process METACOMPASS {
    tag "${meta.id}"
    label 'process_high'

    conda "scgs::metacompass=1.12=1.12--h9948957_4 bioconda::seqkit=2.10.0 bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0e7fe6bd3265990ffcdf96496fe08dc5aa55fd24:62c3137bd1d05677122f8069cb3981ac4e60651e-0' :
        'scgs/mulled-v2-0e7fe6bd3265990ffcdf96496fe08dc5aa55fd24:62c3137bd1d05677122f8069cb3981ac4e60651e-0' }"

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
    # only use reference based
    awk -F'\t' 'NR>1{if(\$3) print \$1}' ${prefix}.metacompass_out/metacompass.tsv > refs.id
    seqkit grep -f refs.id -o ${prefix}.metacompass.ctg.fa ${prefix}.metacompass_out/metacompass.final.ctg.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaCompass: 1.12
    END_VERSIONS
    """
}
