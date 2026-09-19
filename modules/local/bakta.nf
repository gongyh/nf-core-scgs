nextflow.enable.types = true

process BAKTA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bakta=1.11.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.11.4--pyhdfd78af_0' :
        'biocontainers/bakta:1.11.4--pyhdfd78af_0' }"

    input:
    tuple(meta: Map, fasta: Path)
    db: Path
    proteins: List<Path>
    prodigal_tf: List<Path>

    output:
    record(meta: meta, bakta_dir: file("${task.ext.prefix ?: meta.id}", type: "dir"), versions: file("versions.yml"))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_tf_opt = prodigal_tf ? "--prodigal-tf ${prodigal_tf[0]}" : ""
    """
    bakta \\
        $fasta \\
        $args \\
        --threads $task.cpus \\
        --output $prefix \\
        --prefix $prefix \\
        $proteins_opt \\
        $prodigal_tf_opt \\
        --db $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}/${prefix}.embl
    touch ${prefix}/${prefix}.faa
    touch ${prefix}/${prefix}.ffn
    touch ${prefix}/${prefix}.fna
    touch ${prefix}/${prefix}.gbff
    touch ${prefix}/${prefix}.gff3
    touch ${prefix}/${prefix}.hypotheticals.tsv
    touch ${prefix}/${prefix}.hypotheticals.faa
    touch ${prefix}/${prefix}.tsv
    touch ${prefix}/${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """
}
