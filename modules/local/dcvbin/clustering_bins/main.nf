process DCVBIN_BIN {
    tag "${fasta_file.baseName}"

    conda "${moduleDir}/dnaberts.yaml"
    container 'community.wave.seqera.io/library/dnaberts:7a7299083f265248'

    input:
    path vae_features_file
    path cluster_value_file
    path fasta_file

    output:
    path "${fasta_file.baseName}_bins",        emit: bins_dir
    path "${fasta_file.baseName}_prinum.txt",  emit: label_file
    path "${fasta_file.baseName}_scaffolds2bin.tsv", emit: scaffolds2bin
    path "${fasta_file.baseName}_mqc.tsv", emit: mqc_tsv
    path "versions.yml", emit: versions
    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${fasta_file.baseName}"
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
    echo -e "Metric\tValue" > ${prefix}_mqc.tsv
    echo -e "Number of bins\t\${N_BINS}" >> ${prefix}_mqc.tsv
    echo -e "Number of contigs in bins\t\${N_CONTIGS}" >> ${prefix}_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dcvbin: \$(python -c "import dcvbin; print(dcvbin.__version__)" 2>/dev/null || echo "unknown")
    END_VERSIONS
    exit 0
    """
}
