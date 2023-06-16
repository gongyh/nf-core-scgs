process ANEUFINDER {
    label 'process_medium'

    input:
    path("bams/*")
    path("bams/*")

    output:
    path('CNV_output'),  emit: cnv
    path "versions.yml", emit: versions

    when:
    !params.bulk && params.cnv && !single_end && !params.nanopore

    script:
    """
    aneuf.R ./bams CNV_output ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aneufinder: \$(Rscript -e 'v=format(packageVersion("AneuFinder"));cat(v)')
    END_VERSIONS
    """
}
