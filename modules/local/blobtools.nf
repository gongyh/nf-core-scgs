process BLOBTOOLS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::blobtools=1.0.1--py27_3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/blobtools:1.0.1--py27_3' :
        'biocontainers/blobtools:1.0.1--py27_3' }"

    input:
    tuple val(meta), path(contigs)
    tuple val(meta), path(anno)
    tuple val(meta), path(uniprot_anno)
    val has_uniprot
    path db

    output:
    tuple val(meta), path("${prefix}/${prefix}.blobDB*table.txt"), emit: tax
    tuple val(meta), path("${contigs}") ,                           emit: contigs
    tuple val(meta), path("${prefix}") ,                            emit: tax_split
    path "versions.yml",                                            emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def uniprot_anno_cmd = has_uniprot ? "-t $uniprot_anno" : ""
    """
    mkdir -p ${prefix}
    blobtools create -i $contigs -y spades -t $anno $uniprot_anno_cmd -o ${prefix}/${prefix} --db $db
    blobtools view -i ${prefix}/${prefix}.blobDB.json -r all -o ${prefix}/
    blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r phylum --format pdf -o ${prefix}/
    blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r order --format pdf -o ${prefix}/
    blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r family --format pdf -o ${prefix}/
    blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r genus --format pdf -o ${prefix}/
    blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r species --format pdf -o ${prefix}/
    blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r phylum --format png -o ${prefix}/
    blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r order --format png -o ${prefix}/
    blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r family --format png -o ${prefix}/
    blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r genus --format png -o ${prefix}/
    blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r species --format png -o ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blobtools: \$(echo \$(blobtools -v 2>&1) | sed 's/^.*blobtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
