nextflow.enable.types = true

process EGGNOG {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::eggnog-mapper=2.1.11=pyhdfd78af_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.11--pyhdfd78af_0' :
        'biocontainers/eggnog-mapper:2.1.11--pyhdfd78af_0' }"

    input:
    tuple(meta: Map, faa: Path)
    db: Path

    output:
    record(meta: meta, annotations: file("*.emapper.annotations"), mqc_tsv: file("${meta.id}_eggnog_mqc.tsv"))
    topic:
    file("versions.yml") >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    emapper.py -i $faa -o $prefix --data_dir $db --dmnd_db $db/eggnog_proteins.dmnd -m diamond --cpu ${task.cpus}

    cat > "${meta.id}_eggnog_mqc.tsv" <<'EOF'
# id: eggnog
# section_name: EggNOG Annotations
# plot_type: table
EOF
    printf 'Bin\\tAnnotated genes\\tGenes with COG category\\tGenes with KEGG KO\\tGenes with GO terms\\n' >> "${meta.id}_eggnog_mqc.tsv"
    awk -F '\\t' -v bin='${meta.id}' '
        \$1 == "#query" {
            for (i = 1; i <= NF; i++) {
                if (\$i == "COG_category") cog_col = i
                if (\$i == "KEGG_ko") ko_col = i
                if (\$i == "GOs") go_col = i
            }
            next
        }
        !/^#/ && NF > 1 {
            genes++
            if (cog_col && \$(cog_col) != "-" && \$(cog_col) != "") cog++
            if (ko_col && \$(ko_col) != "-" && \$(ko_col) != "") ko++
            if (go_col && \$(go_col) != "-" && \$(go_col) != "") go++
        }
        END { printf "%s\\t%d\\t%d\\t%d\\t%d\\n", bin, genes, cog, ko, go }
    ' ${prefix}.emapper.annotations >> "${meta.id}_eggnog_mqc.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$(echo \$(emapper.py --version | grep emapper 2>&1 ) | cut -d'/' -f1 | sed 's/^.*emapper-//; s/Using.*\$//')
    END_VERSIONS
    """
}
