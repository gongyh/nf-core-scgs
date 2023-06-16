process VG_CONSTRUCT {
    tag "VG_CONSTRUCT"
    label 'process_medium'

    conda "bioconda::vg=1.45.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.45.0--h9ee0642_0':
        'biocontainers/vg:1.45.0--h9ee0642_0' }"

    input:
    path(fasta)
    path(vcf)

    output:
    path("graph.vg"), emit: vg

    script:
    """
    tabix ${vcf}
    vg construct -t ${task.cpus} -r ${fasta} -v ${vcf} > graph.vg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(vg 2>&1 | head -n 1 | sed 's/vg: variation graph tool, version v//;s/ ".*"//' ))
    END_VERSIONS
    """
}
