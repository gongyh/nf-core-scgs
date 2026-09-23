nextflow.enable.types = true

process MMSEQS_TAXONOMY_MULTIQC {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::gawk=5.3.1"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.1' :
        'quay.io/biocontainers/gawk:5.3.1'}"

    input:
    tuple(meta: Map, taxonomy: Path)

    output:
    record(meta: meta, mqc_tsv: file("${meta.id}_mmseqs_taxonomy_mqc.tsv"))

    script:
    """
    cat > "${meta.id}_mmseqs_taxonomy_mqc.tsv" <<'EOF'
# id: mmseqs_taxonomy
# section_name: MMseqs2 Contig Taxonomy
# plot_type: table
EOF
    printf 'Sample\\tAssigned contigs\\tSpecies\\tGenus\\tFamily\\tOrder\\tClass\\tPhylum\\tOther ranks\\n' >> "${meta.id}_mmseqs_taxonomy_mqc.tsv"
    awk -F '\\t' -v sample='${meta.id}' '
        NF >= 4 {
            assigned++
            if (\$3 == "species") species++
            else if (\$3 == "genus") genus++
            else if (\$3 == "family") family++
            else if (\$3 == "order") order++
            else if (\$3 == "class") class++
            else if (\$3 == "phylum") phylum++
            else other++
        }
        END {
            printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\n", sample, assigned, species, genus, family, order, class, phylum, other
        }
    ' ${taxonomy} >> "${meta.id}_mmseqs_taxonomy_mqc.tsv"
    """
}
