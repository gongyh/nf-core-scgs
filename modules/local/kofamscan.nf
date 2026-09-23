nextflow.enable.types = true

process KOFAMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kofamscan=1.3.0 conda-forge::python=3.6.10"
    container "scgs/mulled-v2-ef3cc10895f39bdde312c5e796de361bc231bb29:f6fe8bf9968d952a4cb8cdb90f165e49c40688b8-0"

    input:
    tuple(meta: Map, faa: Path)
    profile: Path
    ko_list: Path

    output:
    record(meta: meta, txt: files("*_KOs_*.txt"), kofamscan: file("*_KOs_ko.kofamscan"), mqc_tsv: file("${meta.id}_kofam_mqc.tsv"))
    topic:
    file("versions.yml") >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} --keep-tabular -o ${prefix}_KOs_detail.txt ${faa}
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} --keep-tabular -r -f mapper -o ${prefix}_KOs_mapper.txt ${faa}
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} --keep-tabular -r -f mapper-one-line -o ${prefix}_KOs_mapper2.txt ${faa}
    kofam_postprocess.py \$(echo \$(which ko_KO.txt)) ${prefix}_KOs_mapper.txt > ${prefix}_KOs_ko.txt
    ln -sf ${prefix}_KOs_ko.txt ${prefix}_KOs_ko.kofamscan
    cat > "${meta.id}_kofam_mqc.tsv" <<'EOF'
# id: kofam
# section_name: KOfam Annotations
# plot_type: table
EOF
    printf 'Bin\\tAnnotated genes\\tKO assignments\\tDistinct KOs\\tDistinct pathways\\n' >> "${meta.id}_kofam_mqc.tsv"
    awk -F '\\t' -v bin='${meta.id}' '
        NF >= 3 {
            assignments++
            genes[\$1] = 1
            split(\$2, pathway_fields, " ")
            pathways[pathway_fields[1]] = 1
            split(\$3, ko_fields, " ")
            kos[ko_fields[1]] = 1
        }
        END {
            for (gene in genes) n_genes++
            for (ko_id in kos) n_kos++
            for (pathway_id in pathways) n_pathways++
            printf "%s\\t%d\\t%d\\t%d\\t%d\\n", bin, n_genes, assignments, n_kos, n_pathways
        }
    ' ${prefix}_KOs_ko.txt >> "${meta.id}_kofam_mqc.tsv"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kofamscan: \$(echo \$(exec_annotation -v 2>&1) | sed 's/^.*exec_annotation //; s/Using.*\$//')
    END_VERSIONS
    """
}
