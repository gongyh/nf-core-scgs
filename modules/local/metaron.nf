nextflow.enable.types = true

process METARON {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    tuple(meta: Map, contigs: Path, gene_model: Path)

    output:
    record(meta: meta, out_operon: file("${prefix}", type: 'dir'))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    metaron.py -n ${prefix} -p op -i ${gene_model} -j ${contigs} -t 2 -o ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaRon: 1.0
    END_VERSIONS
    """
}
