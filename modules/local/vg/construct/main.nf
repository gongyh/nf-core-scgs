process VG_CONSTRUCT {
    label 'process_medium'
    publishDir "${params.outdir}/vg_graph", mode: 'copy'

    conda "bioconda::vg=1.45.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.45.0--h9ee0642_0':
        'biocontainers/vg:1.45.0--h9ee0642_0' }"

    input:
    path fasta
    path graph_vcf

    output:
    path "graph.vg",                              emit: vg

    when:
    params.vcf

    script:
    """
    tabix ${graph_vcf}
    vg construct -t ${task.cpus} -r ${fasta} -v ${graph_vcf} > graph.vg
    """
}
