process FASTQC {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    conda "bioconda::fastqc=0.11.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'biocontainers/fastqc:0.11.9--0' }"

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path("${name}*_fastqc.{zip,html}"),                   emit: result

    script:
    """
    fastqc --threads ${task.cpus} -q $reads
    """
}
