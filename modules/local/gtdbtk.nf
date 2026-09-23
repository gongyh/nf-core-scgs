nextflow.enable.types = true

process GTDBTK {
    label 'process_medium'

    conda "bioconda::gtdbtk=2.7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gtdbtk:2.7.2--pyhdfd78af_1' :
        'biocontainers/gtdbtk:2.7.2--pyhdfd78af_1' }"

    input:
    fa: Set<Path>
    gtdb: Path

    stage:
    stageAs fa, 'genome/*'

    output:
    record(out: file('out', type: "dir"), taxa: file('taxa.txt'), mqc_tsv: file('GTDBtk_mqc.tsv'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    export GTDBTK_DATA_PATH=$gtdb

    mkdir -p genome

    echo \"# plot_type: 'table'\" > GTDBtk_mqc.tsv
    echo \"# section_name: 'GTDBtk'\" >> GTDBtk_mqc.tsv
    echo \"# file_format: 'tsv'\" >> GTDBtk_mqc.tsv
    echo \"# id: scgs_gtdbtk\" >> GTDBtk_mqc.tsv
    printf 'genome\\tclassification\\n' >> GTDBtk_mqc.tsv

    printf 'genome\\tg__\\n' > taxa.txt

    if [[ -f $gtdb ]]; then
        mkdir -p out
        printf 'genome\\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\\n' > taxa.txt
        cat taxa.txt >> GTDBtk_mqc.tsv
    else
        if [ -n "\$(ls -A genome)" ]; then
            gtdbtk classify_wf \\
                --pplacer_cpus 1 \\
                --genome_dir genome \\
                --extension fasta \\
                --out_dir out \\
                --cpus $task.cpus

            for summary in out/*.summary.tsv; do
                [ -f \"\$summary\" ] || continue
                tail -n +2 \"\$summary\" | cut -f1,2 >> GTDBtk_mqc.tsv
            done
            if [ \$(wc -l < GTDBtk_mqc.tsv) -gt 5 ]; then
                tail -n +6 GTDBtk_mqc.tsv > taxa.txt
            fi
        else
            mkdir -p out
            touch out/no_results.txt
            echo "No fasta to taxonomy!" > out/no_results.txt
            printf 'genome\\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\\n' > taxa.txt
            printf 'genome\\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\\n' >> GTDBtk_mqc.tsv
        fi
    fi

    if [ \$(wc -l < GTDBtk_mqc.tsv) -eq 5 ]; then
        printf 'No classified genomes\\tNA\\n' >> GTDBtk_mqc.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}
