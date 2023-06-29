process AUGUSTUS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::augustus=3.5.0=pl5321h700735d_3 tantan=40"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-25b0c981ecfd8d3b08ff5d0fe770fa0aed57e827:2f3083f6f040a1f2ba35c3999b612686446fc7f3-0' :
        'scgs/mulled-v2-25b0c981ecfd8d3b08ff5d0fe770fa0aed57e827:2f3083f6f040a1f2ba35c3999b612686446fc7f3-0' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${prefix}.aa"), emit: faa
    path("${prefix}*")
    path "versions.yml",                   emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # clean id
    cat $contigs | sed 's/_length.*\$//g' > ${prefix}_clean.fasta
    # mask genome
    tantan ${prefix}_clean.fasta > ${prefix}_mask.fasta
    # gene prediction
    augustus --species=${params.augustus_species} --gff3=on --uniqueGeneId=true --protein=on --codingseq=on ${prefix}_mask.fasta > ${prefix}.gff
    # generate proteins
    getAnnoFasta.pl ${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        augustus: 3.5.0
        tantan: 40
    END_VERSIONS
    """
}
