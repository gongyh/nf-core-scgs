process MONOVAR {
    publishDir path: "${pp_outdir}", mode: 'copy',
                saveAs: { filename ->
                    if (filename.indexOf(".vcf") > 0) "$filename" else null }

    input:
    path("*")
    path("*")
    path fa

    output:
    path('monovar.vcf'),  emit: vcf

    when:
    !params.bulk && params.snv && !params.nanopore

    script:
    pp_outdir = "${params.outdir}/monovar"
    """
    ls *.bam > bams.txt
    samtools mpileup -B -d 10000 -q 40 -f $fa -b bams.txt | /opt/MonoVar/src/monovar.py -f $fa -o monovar.vcf -m ${task.cpus} -b bams.txt
    """
}
