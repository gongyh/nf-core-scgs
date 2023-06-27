process STARAMR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda:staramr=0.8.0--pyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staramr:0.8.0--pyhdfd78af_1' :
        'biocontainers/staramr:0.8.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs)
    val(acquired)
    val(point)
    val(species)
    val(only_known)

    output:
    path("${prefix}/*")

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def species = species
    def known_snp = only_known ? "" : "-l 0.4 -r all -u"
    if (acquired && !point) {
    """
    staramr search -o $prefix $contigs
    """
    } else if(point && !acquired) {
    """
    staramr search --pointfinder-organism $species -o $prefix $contigs
    """
    } else {
    """
    staramr search -o $prefix $contigs
    staramr search --pointfinder-organism $species -o $prefix $contigs
    """
    }
}
