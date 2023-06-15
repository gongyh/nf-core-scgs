process SPLIT_CHECKM {
    label 'process_medium'

    input:
    path("results/spades/*")
    path("results/blob/*")
    path("results/prokka/*")
    path("results/kofam/*")
    val split_bac_level
    val split_euk_level

    output:
    path("split/*")

    script:
    """
    cli.py tools scgs_split --level-bacteria ${split_bac_level} --level-eukaryota ${split_euk_level}
    cd split
    samples=(`ls -d *_${split_bac_level}_Bacteria | sed 's/_${split_bac_level}_Bacteria//g'`)
    for sample in \${samples[*]}; do
    mkdir -p \${sample}_${split_bac_level}_checkM
    checkm lineage_wf -t ${task.cpus} -f \${sample}_${split_bac_level}_checkM.txt -x fasta \${sample}_${split_bac_level}_Bacteria \${sample}_${split_bac_level}_checkM || echo "Ignore internal errors!"
    done
    """
}
