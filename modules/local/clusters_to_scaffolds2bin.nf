process CLUSTERS_TO_SCAFFOLDS2BIN {
    tag "clusters_to_scaffolds2bin.nf"
    label 'process_medium'

    conda "conda-forge::coreutils"
    container'community.wave.seqera.io/library/coreutils:9.5--ae99c88a9b28c264'
    input:
    path cluster_file

    output:
    path "scaffolds2bin.tsv", emit: scaffolds2bin
    path "taxvamb_mqc.tsv", emit: mqc_tsv
    path "versions.yml", emit: versions

    script:
    """
    if [ ! -s ${cluster_file} ]; then
        echo "ERROR: cluster file ${cluster_file} is empty or missing" >&2
        exit 1
    fi
    awk -F'\\t' 'NR>1 {print \$2"\t"\$1}' ${cluster_file} > scaffolds2bin.tsv
    N_BINS=\$(cut -f2 scaffolds2bin.tsv | sort -u | wc -l)
    printf "Metric\tValue\\n" > taxvamb_mqc.tsv
    printf "Number of bins recovered\t\${N_BINS}\\n" >> taxvamb_mqc.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clusters_to_scaffolds2bin.nf: \$(echo 1.0)
    END_VERSIONS
    """
}
