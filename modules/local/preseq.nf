process PRESEQ {
    tag "${meta.id}"
    label 'process_single'
    publishDir "${pp_outdir}", mode: 'copy',
                saveAs: { filename ->
                    if (filename.indexOf(".txt") > 0) filename
                    else if (filename.indexOf(".pdf") > 0) filename
                    else null }

    input:
    tuple val(meta), path(sbed)

    output:
    tuple val(meta), path('*.txt'), emit: preseq_for_multiqc
    tuple val(meta), path('*.pdf')

    when:
    !params.nanopore

    script:
    pp_outdir = "${params.outdir}/preseq"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = meta.single_end ? "" : "-P"
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