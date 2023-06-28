process OUTPUT_DOCUMENTATION {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    conda "conda-forge::markdown=3.4.3 conda-forge::pymdown-extensions=10.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9d4085f2843801e3a749ddf5aafb2163e650905b:957aa01b06e937103f54e0d7f72e2ab0c8be9b6f-0' :
        'scgs/mulled-v2-9d4085f2843801e3a749ddf5aafb2163e650905b:957aa01b06e937103f54e0d7f72e2ab0c8be9b6f-0' }"

    input:
    path output_docs

    output:
    path("results_description.html")
    path "versions.yml", emit: versions

    script:
    """
    markdown_to_html.py -o results_description.html $output_docs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        markdown: '3.4.3'
    END_VERSIONS
    """
}
