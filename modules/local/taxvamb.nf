process VAMB_BIN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::vamb=5.0.4"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vamb:5.0.4--pyhdfd78af_0':
        'quay.io/biocontainers/vamb:5.0.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly), path(abundance_tsv), path(bams, stageAs: "bams/*"), path(taxonomy)

    output:
    tuple val(meta), path("${prefix}/scaffolds2bin.tsv")         , emit: scaffolds2bin
    tuple val(meta), path("${prefix}/bins/*.fna.gz")             , emit: bins             , optional: true
    tuple val(meta), path("${prefix}/vae*_clusters_metadata.tsv"), emit: clusters_metadata
    tuple val(meta), path("${prefix}/vae*_clusters_split.tsv")   , emit: clusters_split   , optional: true
    tuple val(meta), path("${prefix}/vae*_clusters_unsplit.tsv") , emit: clusters_unsplit
    tuple val(meta), path("${prefix}/results_taxometer.tsv")     , emit: taxometer_results, optional: true
    tuple val(meta), path("${prefix}/latent.npz")                , emit: latent_encoding  , optional: true
    tuple val(meta), path("${prefix}/abundance.npz")             , emit: abundance
    tuple val(meta), path("${prefix}/composition.npz")           , emit: composition
    tuple val(meta), path("${prefix}/log.txt")                   , emit: log
    tuple val("${task.process}"), val('vamb'), eval("vamb --version | sed 's/Vamb //'"), emit: versions_vamb, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    def mode    = taxonomy ? "taxvamb" : "default"
    depth_input = abundance_tsv ? "--abundance_tsv ${abundance_tsv}" : (bams ? "--bamdir bams/" : error("Neither abundance_tsv nor bams provided"))
    tax_input   = taxonomy ? "--taxonomy ${taxonomy}" : ""
    def min_len = task.ext.min_contig_len ?: 250
    """
    awk -v min=${min_len} 'BEGIN {RS=">"; ORS=""} NR>1 {seq=\$0; gsub(/\\n/, "", seq); if(length(seq) >= min) print ">"\$0}' ${assembly} > filtered.contigs.fasta

    if [ ! -s filtered.contigs.fasta ]; then
        echo "ERROR: No contigs with length >= ${min_len}" >&2
        exit 1
    fi
    grep '^>' filtered.contigs.fasta | sed 's/^>//' > keep_ids.txt
    awk -F'\\t' '
        NR==FNR {
            a[\$1]=\$0;
            next
        }
        FNR==1 {
            print "contigname\\tabundance"
        }
        {
            if (\$1 in a) print a[\$1]
        }
    ' ${abundance_tsv} keep_ids.txt > filtered.abundance.tsv
    vamb bin \\
        ${mode} \\
        -p ${task.cpus} \\
        --outdir ${prefix}/ \\
        --fasta filtered.contigs.fasta \\
        --abundance_tsv filtered.abundance.tsv \\
        ${tax_input} \\
        ${args}

    awk -F'\\t' 'NR>1 {print \$2"\t"\$1}' ${prefix}/vae*_clusters_unsplit.tsv > ${prefix}/scaffolds2bin.tsv
    """

    stub:
    if(bams && abundance_tsv) {
        error("ERROR: Both bams and abundance TSV supplied to Vamb! Please only supply one.")
    }
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/bins

    echo "" | gzip > ${prefix}/bins/${prefix}.1.fna.gz
    echo "" | gzip > ${prefix}/bins/${prefix}.2.fna.gz

    touch ${prefix}/results_taxometer.tsv
    touch ${prefix}/predictor_model.pt
    touch ${prefix}/vae_clusters_metadata.tsv
    touch ${prefix}/vae_clusters_split.tsv
    touch ${prefix}/vae_clusters_unsplit.tsv
    touch ${prefix}/latent.npz
    touch ${prefix}/model.pt
    touch ${prefix}/abundance.npz
    touch ${prefix}/composition.npz
    touch ${prefix}/log.txt
    """
}
