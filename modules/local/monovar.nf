process MONOVAR {
    label 'process_medium'

    input:
    path("*")
    path("*")
    path fa

    output:
    path('monovar.vcf'),  emit: vcf
    path "versions.yml",  emit: versions

    when:
    !params.bulk && params.snv && !params.nanopore

    script:
    pp_outdir = "${params.outdir}/monovar"
    """
    ls *.bam > bams.txt
    samtools mpileup -B -d 10000 -q 40 -f $fa -b bams.txt | /opt/MonoVar/src/monovar.py -f $fa -o monovar.vcf -m ${task.cpus} -b bams.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        monovar: 'v0.0.1'
    END_VERSIONS
    """
}
