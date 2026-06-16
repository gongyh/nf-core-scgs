process CHECKM2 {
    label 'process_medium'

    conda "bioconda::checkm2=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0' :
        'community.wave.seqera.io/library/checkm2:1.0.1--034a3a15afae63b1' }"

    input:
    path input_dir
    val ext
    path db

    output:
    path('checkm2_results.txt'), emit: txt
    path('CheckM2_mqc.tsv')    , emit: mqc_tsv
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def out_file = (ext == "fasta") ? "spades_checkM2.txt" : "checkm2_results.txt"
    """
    checkm2 predict --threads ${task.cpus} -x ${ext} --input ${input_dir} --output-directory checkm2 --database_path ${db}
    cp checkm2/quality_report.tsv ${out_file}
    echo "# plot_type: 'table'" > CheckM2_mqc.tsv
    echo "# section_name: 'CheckM2'" >> CheckM2_mqc.tsv
    cat ${out_file} >> CheckM2_mqc.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm2: \$( checkm2 --version )
    END_VERSIONS
    """
}
