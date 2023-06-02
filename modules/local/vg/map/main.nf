process VG_MAP {
    tag "${prefix}"
    label 'process_low'
    publishDir "${params.outdir}/vg_map", mode: 'copy'

    conda "bioconda::vg=1.45.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.45.0--h9ee0642_0':
        'biocontainers/vg:1.45.0--h9ee0642_0' }"

    input:
    path vg
    path reads
    val(single_end)

    output:
    // file "graph.xg" into graph_xg
    path "*.gam",                            emit: gam
    path "*.stats.txt"

    script:
    prefix = reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    def R1 = reads[0].toString()
    if (single_end) {
    """
    vg index -x graph.xg -g graph.gcsa -k 16 ${vg}
    vg map -t ${task.cpus} -w 1024 -f $R1 -x graph.xg -g graph.gcsa > ${prefix}.gam
    vg stats -z -a ${prefix}.gam ${vg} > stats.txt
    """
    } else {
    def R2 = reads[1].toString()
    """
    vg index -t ${task.cpus} -x graph.xg -g graph.gcsa -k 16 ${vg}
    vg map -t ${task.cpus} -w 1024 -f $R1 $R2 -x graph.xg -g graph.gcsa > ${prefix}.gam
    vg stats -z -a ${prefix}.gam ${vg} > ${prefix}.stats.txt
    """
    }
}
