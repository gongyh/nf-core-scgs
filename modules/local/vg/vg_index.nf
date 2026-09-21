nextflow.enable.types = true

process VG_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vg=1.45.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.45.0--h9ee0642_0':
        'biocontainers/vg:1.45.0--h9ee0642_0' }"

    input:
    tuple(meta: Map, reads: List<Path>, vg: Path)

    output:
    record(meta: meta, gam: file("${prefix}.gam"), txt: file('*.stats.txt'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
    """
    vg index -x graph.xg -g graph.gcsa -k 16 ${vg}
    vg map -t ${task.cpus} -w 1024 -f ${reads[0]} -x graph.xg -g graph.gcsa > ${prefix}.gam
    vg stats -z -a ${prefix}.gam ${vg} > stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(vg 2>&1 | head -n 1 | sed 's/vg: variation graph tool, version v//;s/ ".*"//' ))
    END_VERSIONS
    """
    } else {
    """
    vg index -t ${task.cpus} -x graph.xg -g graph.gcsa -k 16 ${vg}
    vg map -t ${task.cpus} -w 1024 -f ${reads[0]} ${reads[1]} -x graph.xg -g graph.gcsa > ${prefix}.gam
    vg stats -z -a ${prefix}.gam ${vg} > ${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(vg 2>&1 | head -n 1 | sed 's/vg: variation graph tool, version v//;s/ ".*"//' ))
    END_VERSIONS
    """
    }
}
