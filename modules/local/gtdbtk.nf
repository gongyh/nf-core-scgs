process GTDBTK {
    label 'process_medium'

    conda "bioconda::gtdbtk=2.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.1.1--pyhdfd78af_1' :
        'biocontainers/gtdbtk:2.1.1--pyhdfd78af_1' }"

    input:
    path(fa)
    path(gtdb)

    output:
    path("out/*")
    path("genome/*")       , emit: scaffolds
    path('taxa.txt')       , emit: taxa
    path('GTDBtk_mqc.tsv') , emit: mqc_tsv
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export GTDBTK_DATA_PATH=$gtdb

    mkdir -p genome
    cp $fa genome

    echo \"# plot_type: 'table'\" > GTDBtk_mqc.tsv
    echo \"# section_name: 'GTDBtk'\" >> GTDBtk_mqc.tsv

    echo \$'genome\\tg__' >  taxa.txt

    if [[ -f $gtdb ]]; then
        mkdir -p out
        echo \$'genome\\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli' > taxa.txt
    else
        if [ ! -f genome/no_fasta.txt ];then
            gtdbtk classify_wf \\
                --pplacer_cpus 1 \\
                --genome_dir genome \\
                --extension fasta \\
                --out_dir out \\
                --cpus $task.cpus

            if [ -f out/*.summary.tsv ]; then
                cut -f1,2 out/*.summary.tsv | grep -v classification > taxa.txt
                cut -f1,2 out/*.summary.tsv >> GTDBtk_mqc.tsv
            fi
        else
            mkdir -p out
            touch out/no_results.txt
            echo "No fasta to taxonomy!" > out/no_results.txt
            echo \$'genome\\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli' > taxa.txt
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
