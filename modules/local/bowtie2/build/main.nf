process BOWTIE2_BUILD {
    label 'process_high'

    conda "bioconda::bowtie2=2.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.4--py39hbb4e92a_0' :
        'biocontainers/bowtie2:2.4.4--py39hbb4e92a_0' }"

    input:
    path(fasta)

    output:
    path("Bowtie2Index") , emit: bowtie2_index


    when:
    params.fasta

    script:
    """
    mkdir -p Bowtie2Index; cd Bowtie2Index
    ln -s ../${fasta} genome.fa
    bowtie2-build genome.fa genome
    """
}
