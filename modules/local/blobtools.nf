process BLOBTOOLS {
    tag "${prefix}"
    publishDir "${params.outdir}/blob", mode: 'copy'

    input:
    path contigs
    path anno
    val has_uniprot
    path uniprot_anno

    output:
    path("${prefix}/${prefix}.blobDB*table.txt"),              emit: tax
    path("${contigs}"),                                        emit: contigs
    path("${prefix}"),                                         emit: tax_split

    script:
    prefix = contigs.toString() - ~/(\.ctg200\.fasta)?(\.ctg200)?(\.fasta)?(\.fa)?$/
    def uniprot_anno_cmd = has_uniprot ? "-t $uniprot_anno" : ""
    """
    mkdir -p ${prefix}
    blobtools create -i $contigs -y spades -t $anno $uniprot_anno_cmd -o ${prefix}/${prefix} \
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
