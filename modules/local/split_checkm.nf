process SPLIT_CHECKM {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path("results/spades/*")
    path("results/blob/*")
    path("results/prokka/*")
    path("results/kofam/*")

    output:
    path("split/*")

    script:
    def split_bac_level = params.split_bac_level ? params.split_bac_level : "genus"
    def split_euk_level = params.split_euk_level ? params.split_euk_level : "genus"
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
