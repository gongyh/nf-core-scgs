process SPLIT_CHECKM_EUKCC {
    label 'process_medium'

    input:
    path("results/spades/*")
    path("results/blob/*")
    path("results/prokka/*")
    path("results/kofam/*")
    path db
    val split_bac_level
    val split_euk_level

    output:
    path("split/*")

    script:
    """
    export HOME=/tmp/
    if [ -f "/tmp/.etetoolkit/taxa.sqlite" ]; then
    echo "NCBI taxa database exist!"
    else
    python -c "from ete3 import NCBITaxa; ncbi = NCBITaxa(taxdump_file='/opt/nf-core-scgs/taxdump.tar.gz')"
    fi
    cli.py tools scgs_split --level-bacteria ${split_bac_level} --level-eukaryota ${split_euk_level}
    cd split
    samples=(`ls -d *_${split_bac_level}_Bacteria | sed 's/_${split_bac_level}_Bacteria//g'`)
    for sample in \${samples[*]}; do
    mkdir -p \${sample}_${split_bac_level}_checkM
    checkm lineage_wf -t ${task.cpus} -f \${sample}_${split_bac_level}_checkM.txt -x fasta \${sample}_${split_bac_level}_Bacteria \${sample}_${split_bac_level}_checkM || echo "Ignore internal errors!"
    cd \${sample}_${split_euk_level}_Eukaryota
    contigs=(`ls *.fasta`)
    for contig in \${contigs[*]};do
        prefix=\${contig%.fasta}
        # clean id
        cat \$contig | sed 's/_length.*\$//g' > \${prefix}_clean.fasta
        # mask genome
        tantan \${prefix}_clean.fasta > \${prefix}_mask.fasta
        # gene prediction
        augustus --species=${params.augustus_species} --gff3=on --uniqueGeneId=true --protein=on --codingseq=on \${prefix}_mask.fasta > \${prefix}.gff
        # generate proteins
        getAnnoFasta.pl \${prefix}.gff
    done

    if ls *.aa >/dev/null 2>&1;
    then
        faas=(`ls *.aa`)
        for faa in \${faas[*]};do
            faaPrefix=\${faa%.aa}
            eukcc --db ../../$db --ncores ${task.cpus} --outdir \${faaPrefix} --protein \${faa} || echo "Ignore minor errors of eukcc!"
        done
    fi
    cd ../
    done
    """
}
