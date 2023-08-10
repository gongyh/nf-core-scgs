process CHECKM2 {
    label 'process_medium'

    conda "bioconda::checkm2=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0' :
        'biocontainers/checkm2:1.0.1--pyh7cba7a3_0' }"

    input:
    path('spades/*')
    path db

    output:
    path('spades_checkM2.txt'), emit: txt
    path('CheckM2_mqc.tsv')   , emit: mqc_tsv
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    checkm2 predict --threads ${task.cpus} -x fasta --input spades --output-directory checkm2 --database_path db
    cp checkm2/quality_report.tsv spades_checkM2.txt

    echo \"# plot_type: 'table'\" > CheckM2_mqc.tsv
    echo \"# section_name: 'CheckM'\" >> CheckM2_mqc.tsv
    cat spades_checkM2.txt >> CheckM2_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$( checkm2 --version )
    END_VERSIONS
    """
}
