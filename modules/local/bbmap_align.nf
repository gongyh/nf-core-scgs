nextflow.enable.types = true

process BBMAP_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bbmap=39.01 bioconda::samtools=1.16.1 pigz=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0' :
        'biocontainers/mulled-v2-008daec56b7aaf3f162d7866758142b9f889d690:e8a286b2e789c091bac0a57302cdc78aa0112353-0' }"

    input:
    tuple(meta: Map, fastq: List<Path>)
    ref: Path

    output:
    record(meta: meta, clean_fastq: file("*_removehost*.fq.gz"), log: file("*.log"), versions: file("versions.yml"))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input = meta.single_end ? "in=${fastq}" : "in=${fastq[0]} in2=${fastq[1]}"
    def outu = meta.single_end ? "outu=${prefix}_removehost.fq.gz" : "outu=${prefix}_removehost_R1.fq.gz outu2=${prefix}_removehost_R2.fq.gz"

    // Set the db variable to reflect the three possible types of reference input: 1) directory
    // named 'ref', 2) directory named something else (containg a 'ref' subdir) or 3) a sequence
    // file in fasta format
    if ( ref.isDirectory() ) {
        def ref_path = "${ref}"
        if (ref_path ==~ /(.\/)?ref\/?/) {
            db = ''
        } else {
            db = "path=${ref}"
        }
    } else {
        db = "ref=${ref}"
    }

    """
    bbmap.sh \\
        $db \\
        $input \\
        out=${prefix}.bam \\
        $outu \\
        $args \\
        threads=$task.cpus \\
        -Xmx${task.memory.toGiga()}g \\
        &> ${prefix}.bbmap.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
