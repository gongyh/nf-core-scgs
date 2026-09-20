nextflow.enable.types = true

process MONOVAR {
    label 'process_medium'

    conda "bioconda::python-monovar=0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python-monovar:0.1--hdfd78af_0' :
        'biocontainers/python-monovar:0.1--hdfd78af_0' }"

    input:
    bams: Bag<Path>
    bais: Bag<Path>
    fa: Path

    output:
    record(vcf: file('monovar.vcf'))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    """
    ls *.bam > bams.txt
    samtools mpileup -B -d 10000 -q 40 -f $fa -b bams.txt | monovar.py -f $fa -o monovar.vcf -m ${task.cpus} -b bams.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        monovar: 0.0.1
    END_VERSIONS
    """
}
