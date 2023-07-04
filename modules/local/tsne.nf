process TSNE {
    tag "$meta.id"

    conda "conda-forge::opentsne=1.0.0 conda-forge::h5py=3.9.0 conda-forge::numpy=1.25.0 conda-forge::pandas=2.0.2 bioconda::kpal=2.1.1 bioconda::perl-bioperl=1.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8905087433117c98a93e379c07447431e85bdd71:5402918794aa21f8f7e4b46973655d86142c9ffb-0' :
        'scgs/mulled-v2-8905087433117c98a93e379c07447431e85bdd71:5402918794aa21f8f7e4b46973655d86142c9ffb-0' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${prefix}_tsne.tsv"), emit: tsv
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    faFilterByLen.pl ${contigs} 1000 > ${prefix}.ctg1k.fasta
    if [ -s ${prefix}.ctg1k.fasta ]; then
        kpal count -k 4 -r ${prefix}.ctg1k.fasta ${prefix}.4mer
        kmer_tsne.py ${prefix}.4mer ${prefix}_tsne.tsv ${task.cpus}
    else
        touch ${prefix}_tsne.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kpal: \$( checkm 2>&1 | grep 'kpal' | sed 's/.*kpal version//;s/ .*//' )
    END_VERSIONS
    """
}
