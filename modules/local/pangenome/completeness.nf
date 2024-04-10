process COMPLETENESS {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.10.14 bioconda::odgi=0.8.4 bioconda::seqkit=2.8.0 bioconda::gfapy=1.2.3 bioconda::clustalw=2.1 bioconda::gffread=0.12.7 bioconda::aster=1.16 bioconda::blast=2.15.0 conda-forge::biopython=1.83 conda-forge::toytree=3.0.1 conda-forge::pandas=2.2.1 conda-forge::numpy=1.26.4 conda-forge::prettytable=3.10.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-19582d58b6424cdcb44cb604f72b86d0a8ee2ce1:7ddd327f257df15af9179f8c340d2bc793185cb0-0':
        'scgs/mulled-v2-19582d58b6424cdcb44cb604f72b86d0a8ee2ce1:7ddd327f257df15af9179f8c340d2bc793185cb0-0' }"

    input:
    tuple val(meta), path(contigs)
    val(genusName)
    path(pangenomeDB)
    path(coreGenesFile)

    output:
    path("${prefix}/completeness.txt")  , emit: txt

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mgpgtools.py core -db $pangenomeDB -name $genusName -fasta $contigs -coreGenes $coreGenesFile -outdir ${prefix}
    """
}
