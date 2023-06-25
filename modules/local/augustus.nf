process AUGUSTUS {
    tag "$meta.id"
    label 'process_medium'

    conda "augustus=3.5.0--pl5321h700735d_3,tantan=40,perl-getopt-long=2.54"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-383531f26ddbe208b657a2b4a8d48b7c3644bdc1:d12534b5980a00347a7ae2723726dc1a81248f12-0' :
        'scgs/mulled-v2-383531f26ddbe208b657a2b4a8d48b7c3644bdc1:d12534b5980a00347a7ae2723726dc1a81248f12-0' }"

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${prefix}.aa"), emit: faa
    path("${prefix}*")

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
    """
}
