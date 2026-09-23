nextflow.enable.types = true

process DCVBIN_BIN {
    tag "$meta.id"

    conda "${moduleDir}/dnaberts.yaml"
    container 'community.wave.seqera.io/library/dnaberts:7a7299083f265248'

    input:
    tuple(meta: Map, vae_features_file: Path, cluster_value_file: Path, fasta_file: Path)

    output:
    record(meta: meta, bins_dir: file("${prefix}_bins", type: 'dir'), label_file: file("${prefix}_prinum.txt"), scaffolds2bin: file("${prefix}_scaffolds2bin.tsv"), mqc_tsv: file("${prefix}_mqc.tsv"))
    topic:
    file('versions.yml') >> 'versions'

    script:
    def args    = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    set +e
    export QT_XCB_GL_INTEGRATION="none"
    mkdir -p "${prefix}_bins"
    python ${projectDir}/bin/dcvbin/scripts/myCluster_2.py \
        -vd "${vae_features_file}" \
        -ld "${prefix}_prinum.txt" \
        -fd "${fasta_file}" \
        -bd "${prefix}_bins" \
        -cvf "${cluster_value_file}"
    grep '^>' "${fasta_file}" | sed 's/^>//' > ${prefix}_contig_names.txt
    paste ${prefix}_contig_names.txt ${prefix}_prinum.txt > ${prefix}_scaffolds2bin.tsv
    N_BINS=\$(ls -1 ${prefix}_bins/*.fa 2>/dev/null | wc -l)
    N_CONTIGS=\$(wc -l < ${prefix}_prinum.txt)
    echo -e "Metric\\tValue" > ${prefix}_mqc.tsv
    echo -e "Number of bins\\t\${N_BINS}" >> ${prefix}_mqc.tsv
    echo -e "Number of contigs in bins\\t\${N_CONTIGS}" >> ${prefix}_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dcvbin: 732ee4257f7da200994f4c105e9dacbc74242883
    END_VERSIONS
    exit 0
    """
}
