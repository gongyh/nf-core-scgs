process SPLIT_CHECKM {
    label 'process_medium'

    conda "bioconda::checkm-genome=1.2.1 bioconda::augustus=3.5.0=pl5321h700735d_3 bioconda::tantan=40 conda-forge::biopython=1.81 conda-forge::typer=0.9.0 conda-forge::numpy=1.25.0 bioconda::eukcc=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-28c5d03d1ac8475499ba2a43715feecc3e991223:c795f73b9d282e25900663d2b634c26711c5b8a4-0' :
        'scgs/mulled-v2-28c5d03d1ac8475499ba2a43715feecc3e991223:c795f73b9d282e25900663d2b634c26711c5b8a4-0' }"

    input:
    tuple val(meta), path(contigs)
    tuple val(meta), path(blob, stageAs:"blob_")
    tuple val(meta), path(prokka, stageAs:"prokka_")
    tuple val(meta), path(kofam)
    val split_bac_level
    val split_euk_level

    output:
    path("split/*")          , emit: out_put
    path "split/versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prokka_exist = (prokka == null ? false : true)
    def kofam_exist = (kofam == null ? false : true)
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ ! -d results ];then
        mkdir -p results/spades results/blob
        cp $contigs results/spades
        cp -R blob_ results/blob
        mv results/blob/blob_ results/blob/${prefix}
        if [ ${prokka_exist} ];then
            mkdir -p results/prokka
            cp -R prokka_ results/prokka
            mv results/prokka/prokka_ results/prokka/${prefix}
        fi
        if [ ${kofam_exist} ];then
            mkdir -p results/kofam
            cp $kofam results/kofam
        fi
    fi
    cli-single.py tools scgs_split --level-bacteria ${split_bac_level} --level-eukaryota ${split_euk_level} --sample ${prefix}
    cd split
    mkdir -p \${prefix}_${split_bac_level}_checkM
    checkm lineage_wf -t ${task.cpus} -f \${prefix}_${split_bac_level}_checkM.txt -x fasta \${prefix}_${split_bac_level}_Bacteria \${prefix}_${split_bac_level}_checkM || echo "Ignore internal errors!"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$( checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
    END_VERSIONS
    """
}
