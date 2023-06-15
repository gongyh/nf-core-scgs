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
    tuple val(meta), path('*.txt'), emit: txt
    tuple val(meta), path('*.pdf')
    path "versions.yml",            emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """
    } else {
    """
    preseq c_curve ${mode} -s 1e+5 -o ${prefix}_c.txt $sbed
    preseq lc_extrap ${mode} -s 1e+5 -D -o ${prefix}_lc.txt $sbed
    plotPreSeq.R ${prefix}_lc.txt ${prefix}_lc
    preseq gc_extrap -w 1000 -s 1e+7 -D -o ${prefix}_gc.txt $sbed
    plotPreSeq.R ${prefix}_gc.txt ${prefix}_gc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """
    }
}
