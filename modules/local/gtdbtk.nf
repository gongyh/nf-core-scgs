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
    path("genome/*.fasta"), emit: scaffolds
    path('taxa.txt')      , emit: taxa
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export GTDBTK_DATA_PATH=$gtdb

    mkdir -p genome
    cp $fa genome
    echo 'g__' >  taxa.txt

    if [[ -f $gtdb ]]; then
        mkdir -p out
        echo 'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli' > taxa.txt
    else
        gtdbtk classify_wf \\
            --pplacer_cpus 1 \\
            --genome_dir genome \\
            --extension fasta \\
            --out_dir out \\
            --cpus $task.cpus

        if [ -f out/*.summary.tsv ]; then
            cut -f2 out/*.summary.tsv | grep -v classification > taxa.txt
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
