process KRAKEN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kraken2=2.1.2 bioconda::krona=2.7.1 bioconda::krakentools=1.2 bioconda::bracken=3.1"
    container "scgs/mulled-v2-2e2a18ac791581ea95fced5830f3fe8013145898:c5d1b87c47ed8c1dcf991ed390fb3bf63b5342f8-0"

    input:
    tuple val(meta), path(reads)
    path db
    path taxonomy, stageAs: 'taxonomy.tab'

    output:
    tuple val(meta), path("*.krk")   , emit: report
    tuple val(meta), path("*.html")  , emit: html
    path("${prefix}.TDA_genus.txt")  , emit: tda
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode = meta.single_end ? "" : "--paired"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    TAXONOMY=\$(find -L . -name '*.tab' -exec dirname {} \\;)
    kraken2 --db $db --threads ${task.cpus} --report ${prefix}.krk --output ${prefix}.k2 --gzip-compressed ${mode} $reads
    kreport2krona.py -r ${prefix}.krk -o ${prefix}.krn
    ktImportText -o ${prefix}_taxonomy.html ${prefix}.krn
    # Taxonomic Discovery Algorithm
    bracken -d $db -i ${prefix}.krk -o ${prefix}.bracken -w /dev/null -r 150 -l G
    awk -F '\\t' 'BEGIN{ print "genus\tabundance" }NR>1{print \$1"|"\$2"\\t"\$7}' ${prefix}.bracken > ${prefix}.TDA_genus.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken --version 2>&1) | sed 's/^.*kraken //; s/Using.*\$//')
    END_VERSIONS
    """
}
