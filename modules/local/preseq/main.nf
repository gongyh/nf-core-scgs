process PRESEQ {
    tag "${prefix}"
    label 'process_single'
    publishDir "${pp_outdir}", mode: 'copy',
                saveAs: { filename ->
                    if (filename.indexOf(".txt") > 0) filename
                    else if (filename.indexOf(".pdf") > 0) filename
                    else null }

    conda "bioconda::preseq=3.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/preseq:3.1.2--h445547b_2':
        'biocontainers/preseq:3.1.2--h445547b_2' }"

    input:
    path sbed
    val(single_end)

    output:
    path('*.txt'), emit: preseq_for_multiqc
    path('*.pdf')

    when:
    !params.nanopore

    script:
    pp_outdir = "${params.outdir}/preseq"
    prefix = sbed.toString() - ~/(\.markdup\.bed)?(\.markdup)?(\.bed)?$/
    def mode = single_end ? "" : "-P"
    if (params.bulk) {
    """
    preseq c_curve ${mode} -s 1e+5 -o ${prefix}_c.txt $sbed
    preseq lc_extrap ${mode} -s 1e+5 -D -o ${prefix}_lc.txt $sbed
    plotPreSeq.R ${prefix}_lc.txt ${prefix}_lc
    """
    } else {
    """
    preseq c_curve ${mode} -s 1e+5 -o ${prefix}_c.txt $sbed
    preseq lc_extrap ${mode} -s 1e+5 -D -o ${prefix}_lc.txt $sbed
    plotPreSeq.R ${prefix}_lc.txt ${prefix}_lc
    preseq gc_extrap -w 1000 -s 1e+7 -D -o ${prefix}_gc.txt $sbed
    plotPreSeq.R ${prefix}_gc.txt ${prefix}_gc
    """
    }
}
