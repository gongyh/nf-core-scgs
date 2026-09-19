nextflow.enable.types = true

process KRAKEN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kraken2=2.1.2 bioconda::krona=2.7.1 bioconda::krakentools=1.2 bioconda::bracken=3.1"
    container "scgs/mulled-v2-2e2a18ac791581ea95fced5830f3fe8013145898:c5d1b87c47ed8c1dcf991ed390fb3bf63b5342f8-0"

    input:
    tuple(meta: Map, reads: List<Path>)
    db: Path
    taxonomy: Path

    output:
    record(meta: meta, report: file("*.krk"), html: file("*.html"), tda: file("*.TDA_genus.txt"), versions: file("versions.yml"))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    def mode = meta.single_end ? "" : "--paired"
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_args = reads.join(' ')
    """
    TAXONOMY=\$(find -L . -name '*.tab' -exec dirname {} \\;)
    kraken2 --db $db --threads ${task.cpus} --report ${prefix}.krk --output ${prefix}.k2 --gzip-compressed ${mode} ${read_args}
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
