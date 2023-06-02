process OUTPUT_DOCUMENTATION {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path output_docs

    output:
    path("results_description.html")

    script:
    """
    markdown_to_html.py -o results_description.html $output_docs
    """
}
