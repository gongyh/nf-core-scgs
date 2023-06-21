process POINTFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda:staramr==0.8.0--pyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/staramr:0.8.0--pyhdfd78af_1' :
        'biocontainers/staramr:0.8.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contigs)
    path db

    output:
    path("${prefix}/*")

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def species = params.pointfinder_species
    def known_snp = params.only_known ? "" : "-l 0.4 -r all -u"
    """
    mkdir -p $prefix
    python /opt/pointfinder/PointFinder.py -p $db \
    -m blastn -m_p /opt/conda/bin/blastn $known_snp \
    -i $contigs -o $prefix -s $species
    rm -rf $prefix/tmp

    // staramr search --pointfinder-organism $species -o $prefix $contigs
    """
}
