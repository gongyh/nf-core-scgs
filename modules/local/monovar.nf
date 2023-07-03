process MONOVAR {
    label 'process_medium'

    conda "conda-forge::numpy=1.25.0 conda-forge::scipy=1.10.1 bioconda::pysam=0.21.0 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python-monovar:0.1--hdfd78af_0' :
        'biocontainers/python-monovar:0.1--hdfd78af_0' }"

    input:
    path("*")
    path("*")
    path fa

    output:
    path('monovar.vcf'), emit: vcf
    path "versions.yml", emit: versions

    when:
    !params.bulk && params.snv && !params.nanopore

    script:
    pp_outdir = "${params.outdir}/monovar"
    """
    ls *.bam > bams.txt
    samtools mpileup -B -d 10000 -q 40 -f $fa -b bams.txt | monovar_cli.py -f $fa -o monovar.vcf -m ${task.cpus} -b bams.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        monovar: 0.0.1
    END_VERSIONS
    """
}
