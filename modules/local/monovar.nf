process MONOVAR {
    label 'process_medium'

    conda "bioconda::python-monovar=0.1"
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
    task.ext.when == null || task.ext.when

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
