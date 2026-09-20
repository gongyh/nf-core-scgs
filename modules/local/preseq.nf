nextflow.enable.types = true

process PRESEQ {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::preseq=3.2.0 conda-forge::r-base=4.3.0"
    container "scgs/mulled-v2-f75ca76f6f0d8dac03a420a64d4d702441604c14:03f4a075e359bb32a613b098d13dba7b4c8c967f-0"

    input:
    tuple(meta: Map, sbed: Path)

    output:
    record(meta: meta, results: file('preseq', type: 'dir'))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    pp_outdir = "${params.outdir}/preseq"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = meta.single_end ? "" : "-P"
    if (params.bulk) {
    """
    preseq c_curve ${mode} -s 1e+5 -o ${prefix}_c.txt $sbed
    preseq lc_extrap ${mode} -s 1e+5 -D -o ${prefix}_lc.txt $sbed
    plotPreSeq.R ${prefix}_lc.txt ${prefix}_lc

    mkdir -p preseq
    for result in *.txt *.pdf; do
        [ -e "\$result" ] && mv "\$result" preseq/
    done

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

    mkdir -p preseq
    for result in *.txt *.pdf; do
        [ -e "\$result" ] && mv "\$result" preseq/
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """
    }
}
