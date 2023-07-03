process SPLIT_CHECKM {
    label 'process_medium'

    conda "bioconda::checkm-genome=1.2.1 bioconda::augustus=3.5.0=pl5321h700735d_3 bioconda::tantan=40 conda-forge::biopython=1.81 conda-forge::typer=0.9.0 conda-forge::numpy=1.25.0 bioconda::eukcc=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-28c5d03d1ac8475499ba2a43715feecc3e991223:c795f73b9d282e25900663d2b634c26711c5b8a4-0' :
        'scgs/mulled-v2-28c5d03d1ac8475499ba2a43715feecc3e991223:c795f73b9d282e25900663d2b634c26711c5b8a4-0' }"

    input:
    path("results/spades/*")
    path("results/blob/*")
    path("results/prokka/*")
    path("results/kofam/*")
    val split_bac_level
    val split_euk_level

    output:
    path("split/*")
    path "split/versions.yml", emit: versions

    script:
    """
    cli.py tools scgs_split --level-bacteria ${split_bac_level} --level-eukaryota ${split_euk_level}
    cd split
    samples=(`ls -d *_${split_bac_level}_Bacteria | sed 's/_${split_bac_level}_Bacteria//g'`)
    for sample in \${samples[*]}; do
        mkdir -p \${sample}_${split_bac_level}_checkM
        checkm lineage_wf -t ${task.cpus} -f \${sample}_${split_bac_level}_checkM.txt -x fasta \${sample}_${split_bac_level}_Bacteria \${sample}_${split_bac_level}_checkM || echo "Ignore internal errors!"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$( checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
    END_VERSIONS
    """
}
