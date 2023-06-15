process VG_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vg=1.45.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.45.0--h9ee0642_0':
        'biocontainers/vg:1.45.0--h9ee0642_0' }"

    input:
    path(vg)
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.gam"), emit: gam
    path("*.stats.txt")
    path "versions.yml", emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
    """
    vg index -x graph.xg -g graph.gcsa -k 16 ${vg}
    vg map -t ${task.cpus} -w 1024 -f ${reads[0]} -x graph.xg -g graph.gcsa > ${prefix}.gam
    vg stats -z -a ${prefix}.gam ${vg} > stats.txt
    """
    } else {
    """
    vg index -t ${task.cpus} -x graph.xg -g graph.gcsa -k 16 ${vg}
    vg map -t ${task.cpus} -w 1024 -f ${reads[0]} ${reads[1]} -x graph.xg -g graph.gcsa > ${prefix}.gam
    vg stats -z -a ${prefix}.gam ${vg} > ${prefix}.stats.txt
    """
    }
}
