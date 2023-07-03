process PRESEQ {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::preseq=3.2.0 conda-forge::r-base=4.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f75ca76f6f0d8dac03a420a64d4d702441604c14:03f4a075e359bb32a613b098d13dba7b4c8c967f-0':
        'scgs/mulled-v2-f75ca76f6f0d8dac03a420a64d4d702441604c14:03f4a075e359bb32a613b098d13dba7b4c8c967f-0' }"

    input:
    tuple val(meta), path(sbed)

    output:
    tuple val(meta), path('*.txt'), emit: txt
    tuple val(meta), path('*.pdf')
    path "versions.yml"           , emit: versions

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
    preseq gc_extrap -w 1000 -s 1e+7 -B -D -o ${prefix}_gc.txt $sbed
    plotPreSeq.R ${prefix}_gc.txt ${prefix}_gc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """
    }
}
