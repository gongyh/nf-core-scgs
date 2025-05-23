process SPLIT_CHECKM_EUKCC {
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
    path db
    val split_bac_level
    val split_euk_level

    output:
    path("split/*")            , emit: output
    path("split/fa/*")         , emit: fa
    path("split/*.csv")        , emit: csv
    path("split/versions.yml") , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cli.py tools scgs_split --level-bacteria ${split_bac_level} --level-eukaryota ${split_euk_level}
    cd split
    if [ ! -d fa ];then
        mkdir fa
    fi
    samples=(`ls -d *_${split_bac_level}_Bacteria | sed 's/_${split_bac_level}_Bacteria//g'`)
    for sample in \${samples[*]}; do
        mkdir -p \${sample}_${split_bac_level}_checkM
        checkm lineage_wf -t ${task.cpus} --tab_table -f \${sample}_${split_bac_level}_checkM.tsv -x fasta \${sample}_${split_bac_level}_Bacteria \${sample}_${split_bac_level}_checkM || echo "Ignore internal errors!"
        cd \${sample}_${split_euk_level}_Eukaryota
        if ls *.fasta >/dev/null 2>&1;
        then
            contigs=(`ls *.fasta`)
            for contig in \${contigs[*]};do
                prefix=\${contig%.fasta}
                # clean id
                cat \$contig | sed 's/_length.*\$//g' > \${prefix}_clean.fasta
                eukcc single --db ../../$db --out \${prefix} --threads ${task.cpus} \${prefix}_clean.fasta || echo "Ignore minor errors of eukcc!"
            done
        fi
        cd ../
    done
    prepare_fa.py
    if [ "`ls -A fa`" == "" ];then
        touch fa/no_fasta.txt
        echo "No splited fasta with integrity greater than 40% and contamination less than 10%" > fa/no_fasta.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        augustus: 3.5.0
        tantan: 40
        checkm: \$( checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
        eukcc: \$(echo \$(eukcc -v 2>&1) | sed 's/^.*EukCC version //; s/Using.*\$//')
    END_VERSIONS
    """
}
