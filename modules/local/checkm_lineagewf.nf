process CHECKM_LINEAGEWF {
    label 'process_medium'

    conda "bioconda::checkm-genome=1.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkm-genome:1.2.1--pyhdfd78af_0' :
        'biocontainers/checkm-genome:1.2.1--pyhdfd78af_0' }"

    input:
    path('spades/*')
    val(genus)

    output:
    path('spades_checkM.txt')
    path "versions.yml",       emit: versions

    script:
    def checkm_wf = genus ? "taxonomy_wf" : "lineage_wf"
    """
    if [ \"${checkm_wf}\" == \"taxonomy_wf\" ]; then
    checkm taxonomy_wf -t ${task.cpus} --tab_table -f spades_checkM.txt -x fasta genus ${params.genus} spades spades_checkM
    else
    checkm lineage_wf -t ${task.cpus} -r --tab_table -f spades_checkM.txt -x fasta spades spades_checkM
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$( checkm 2>&1 | grep '...:::' | sed 's/.*CheckM v//;s/ .*//' )
    END_VERSIONS
    """
}
