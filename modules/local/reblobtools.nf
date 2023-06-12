process REBLOBTOOLS {
    tag "$meta.id"

    input:
    tuple val(meta), path contigs
    tuple val(meta), path anno
    val has_uniprot
    tuple val(meta), path uniprot_anno
    path bam

    output:
    tuple val(meta), path("${prefix}/${prefix}.blobDB*table.txt")
    tuple val(meta), path("${contigs}")
    tuple val(meta), path("${prefix}")

    when:
    params.remap

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def uniprot_anno_cmd = has_uniprot ? "-t $uniprot_anno" : ""
    """
    mkdir -p ${prefix}
        samtools sort -o ${prefix}_ass.sort.bam ${prefix}_ass.bam
        samtools index ${prefix}_ass.sort.bam
    blobtools create -i $contigs -y spades -t $anno $uniprot_anno_cmd -b ${prefix}_ass.sort.bam -o ${prefix}/${prefix} \
    --db /opt/conda/envs/nf-core-gongyh-scgs/lib/python3.6/site-packages/data/nodesDB.txt
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
    """
}
