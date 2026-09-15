process BAKTA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bakta=1.11.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bakta:1.11.4--pyhdfd78af_0' :
        'biocontainers/bakta:1.11.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("${prefix}/${prefix}.embl")             , emit: embl
    tuple val(meta), path("${prefix}/${prefix}.faa")              , emit: faa
    tuple val(meta), path("${prefix}/${prefix}.ffn")              , emit: ffn
    tuple val(meta), path("${prefix}/${prefix}.fna")              , emit: fna
    tuple val(meta), path("${prefix}/${prefix}.gbff")             , emit: gbff
    tuple val(meta), path("${prefix}/${prefix}.gff3")             , emit: gff
    tuple val(meta), path("${prefix}/${prefix}.hypotheticals.tsv"), emit: hypotheticals_tsv
    tuple val(meta), path("${prefix}/${prefix}.hypotheticals.faa"), emit: hypotheticals_faa
    tuple val(meta), path("${prefix}/${prefix}.tsv")              , emit: tsv
    tuple val(meta), path("${prefix}/${prefix}.txt")              , emit: txt
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_tf = prodigal_tf ? "--prodigal-tf ${prodigal_tf[0]}" : ""
    """
    bakta \\
        $fasta \\
        $args \\
        --threads $task.cpus \\
        --output $prefix \\
        --prefix $prefix \\
        $proteins_opt \\
        $prodigal_tf \\
        --db $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
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
