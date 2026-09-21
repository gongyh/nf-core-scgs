nextflow.enable.types = true

process VAMB_BIN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::vamb=5.0.4"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vamb:5.0.4--pyhdfd78af_0':
        'quay.io/biocontainers/vamb:5.0.4--pyhdfd78af_0' }"

    input:
    tuple(meta: Map, assembly: Path, abundance_tsv: Path, taxonomy: Path)

    output:
    record(meta: meta, scaffolds2bin: file("${prefix}/scaffolds2bin.tsv"), bins: files("${prefix}/bins/*.fna.gz", optional: true), clusters_metadata: file("${prefix}/vae*_clusters_metadata.tsv"), clusters_split: file("${prefix}/vae*_clusters_split.tsv", optional: true), clusters_unsplit: file("${prefix}/vae*_clusters_unsplit.tsv"), taxometer_results: file("${prefix}/results_taxometer.tsv", optional: true), latent_encoding: file("${prefix}/latent.npz", optional: true), abundance: file("${prefix}/abundance.npz"), composition: file("${prefix}/composition.npz"), log: file("${prefix}/log.txt"))
    topic:
    file('versions.yml') >> 'versions'

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    def mode    = "taxvamb"
    tax_input   = "--taxonomy ${taxonomy}"
    def min_len = task.ext.min_contig_len ?: '250'
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vamb: \$(vamb --version | sed 's/Vamb //')
    END_VERSIONS
    """

    stub:
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
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vamb: stub
    END_VERSIONS
    """
}
