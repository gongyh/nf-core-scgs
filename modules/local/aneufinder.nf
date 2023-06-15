process ANEUFINDER {
    publishDir "${pp_outdir}", mode: 'copy'

    input:
    path("bams/*")
    path("bams/*")

    output:
    path('CNV_output'),  emit: cnv
    path "versions.yml", emit: versions

    when:
    !params.bulk && params.cnv && !single_end && !params.nanopore

    script:
    pp_outdir = "${params.outdir}/aneufinder"
    """
    aneuf.R ./bams CNV_output ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aneufinder: \$(Rscript -e 'print(packageVersion("AneuFinder"))' | sed 's/^.*AneuFinder //;')
    END_VERSIONS
    """
}
