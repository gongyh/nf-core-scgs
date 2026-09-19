nextflow.enable.types = true

process ANEUFINDER {
    label 'process_medium'

    conda "bioconda::bioconductor-aneufinder=1.26.0=r42hf17093f_1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aneufinder:1.26.0--r42hf17093f_1' :
        'biocontainers/bioconductor-aneufinder:1.26.0--r42hf17093f_1' }"

    input:
    bams: Bag<Path>
    bais: Bag<Path>

    output:
    record(cnv: file('CNV_output', type: 'dir'), versions: file('versions.yml'))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    """
    mkdir bams
    cd bams && ln -s ../*.bam ../*.bai . && cd ..
    aneuf.R ./bams CNV_output ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aneufinder: \$(Rscript -e 'v=format(packageVersion("AneuFinder"));cat(v)')
    END_VERSIONS
    """
}
