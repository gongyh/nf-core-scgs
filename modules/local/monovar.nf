process MONOVAR {
    label 'process_medium'

    conda "conda-forge::numpy=1.25.0 conda-forge::scipy=1.10.1 bioconda::pysam=0.21.0 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-3fac00b995a603e53b168c14fd623583081a1b9d:dc88b89d94c822ade5a35acdb5836139bb931890-0' :
        'scgs/mulled-v2-3fac00b995a603e53b168c14fd623583081a1b9d:dc88b89d94c822ade5a35acdb5836139bb931890-0' }"

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
    samtools mpileup -B -d 10000 -q 40 -f $fa -b bams.txt | monovar_cli.py -f $fa -o monovar.vcf -m ${task.cpus} -b bams.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        monovar: v0.0.1
    END_VERSIONS
    """
}
