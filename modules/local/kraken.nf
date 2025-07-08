process KRAKEN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kraken2=2.1.2 bioconda::krona=2.7.1 bioconda::krakentools=1.2"
    container "scgs/mulled-v2-3bbb1b9ff2130265cf8d9498a097b04978fb988f:6688dcb6662e35001e709b425821fff321f15540-0"

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    TAXONOMY=\$(find -L . -name '*.tab' -exec dirname {} \\;)
    kraken2 --db $db --threads ${task.cpus} --report ${prefix}.krk --output ${prefix}.k2 --gzip-compressed ${mode} $reads
    kreport2krona.py -r ${prefix}.krk -o ${prefix}.krn
    ktImportText -o ${prefix}_taxonomy.html ${prefix}.krn
    # Taxonomic Discovery Algorithm
    cat ${prefix}.krn | grep f__ | grep g__ | grep -v s__ > genus_${prefix}.krn
    total_sum=\$(awk -F '\t' '{sum += \$1} END {print sum}' genus_${prefix}.krn)
    awk -F '\t' -v total="\$total_sum" '
    BEGIN{ print "genus\tabundance" }
    {
        f_val = "";
        g_val = "";
        for (i=1; i<=NF; i++) {
            if (\$i ~ /^f__/) {
                f_val = \$i;
            }
            if (\$i ~ /^g__/) {
                g_val = \$i;
            }
        }
        percent = (\$1 / total) * 100;
        printf "%s|%s\t%.2f\n", f_val, g_val, percent;
    }' genus_${prefix}.krn > ${prefix}.TDA_genus.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken --version 2>&1) | sed 's/^.*kraken //; s/Using.*\$//')
    END_VERSIONS
    """
}
