nextflow.enable.types = true

process GET_SOFTWARE_VERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    versions: Path

    output:
    record(yml: file('software_versions.yml'), mqc_yml: file('software_versions_mqc.yml'), versions: file('versions.yml'))

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
