process VG_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vg=1.45.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.45.0--h9ee0642_0':
        'biocontainers/vg:1.45.0--h9ee0642_0' }"

    input:
    path(vg)
    tuple val(meta), path(gam)

    output:
    path("*.calls.vcf"), emit: call

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vg augment -t ${task.cpus} ${vg} ${prefix}.gam -A ${prefix}.aug.gam > ${prefix}.aug.vg
    vg index -t ${task.cpus} ${prefix}.aug.vg -x ${prefix}.aug.xg
    vg pack -t ${task.cpus} -x ${prefix}.aug.xg -g ${prefix}.aug.gam -Q 5 -s 5 -o ${prefix}.aln_aug.pack
    vg call -t ${task.cpus} ${prefix}.aug.xg -k ${prefix}.aln_aug.pack > ${prefix}.calls.vcf
    """
}
