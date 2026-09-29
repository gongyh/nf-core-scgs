nextflow.enable.types = true

process METABULI_TAXONOMY_TABLES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.1"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.1' :
        'quay.io/biocontainers/gawk:5.3.1'}"

    input:
    tuple(meta: Map, classifications: Path)

    output:
    record(meta: meta, taxonomy: file("${meta.id}_taxvamb_taxonomy.tsv"), semibin_taxonomy: file("${meta.id}_semibin_taxonomy.tsv"), mqc_tsv: file("${meta.id}_metabuli_taxonomy_mqc.tsv"))

    script:
    """
    awk -F '\\t' -v taxvamb='${meta.id}_taxvamb_taxonomy.tsv' -v semibin='${meta.id}_semibin_taxonomy.tsv' '
        BEGIN {OFS="\\t"; print "contigs", "predictions" > taxvamb}
        /^#/ {next}
        {
            if (NF < 8) {
                print "ERROR: Metabuli classifications require --lineage 1" > "/dev/stderr"
                exit 1
            }
            lineage = (\$1 == 0 || \$8 == "" || \$8 == "-") ? "unknown" : \$8
            print \$2, lineage > taxvamb
            count = split(lineage, taxa, ";")
            name = taxa[count]
            sub(/^[^_]*_/, "", name)
            normalized = ""
            for (i = 1; i <= count; i++) {
                sub(/^[^_]*_/, "", taxa[i])
                normalized = normalized (i > 1 ? ";" : "") taxa[i]
            }
            rank = (lineage == "unknown") ? "no rank" : \$7
            print \$2, \$3, rank, name, 0, 0, 0, \$5, normalized > semibin
            total++
            if (lineage != "unknown") assigned++
        }
        END {
            print "# id: metabuli_taxonomy"
            print "# section_name: Metabuli Contig Taxonomy"
            print "# plot_type: table"
            print "Sample", "Contigs", "Assigned contigs", "Unclassified contigs"
            print "${meta.id}", total+0, assigned+0, total-assigned
        }
    ' ${classifications} > ${meta.id}_metabuli_taxonomy_mqc.tsv
    touch ${meta.id}_semibin_taxonomy.tsv
    """
}
