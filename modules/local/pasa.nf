process PANTA {
    tag "panta"
    label 'process_high'

    conda "python=3.7.10 biopython=1.79 bedtools=2.30.0 prodigal=2.6.3 cd-hit=4.8.1 blast=2.13.0 hmmer=3.3.2 diamond=2.0.14 mcl=14.137 mafft=7.526 parallel=20220222 numpy=1.21.6 scipy=1.7.3 numba=0.56.3 networkx=2.6.3 pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-adf032e06d3ecc4986e4c4efa160e481707acd5a:5d5dd40d2a8d9df574dfc603ccb059d9d3354e19-0' :
        'scgs/mulled-v2-adf032e06d3ecc4986e4c4efa160e481707acd5a:5d5dd40d2a8d9df574dfc603ccb059d9d3354e19-0' }"

    input:
    path(refs_fna)

    output:
    path("panta_refs")                , emit: db
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p gffs
    refs=(${refs_fna})
    for fna in \${refs[*]}; do
        if [[ \$fna == *.fna.gz ]]; then
            bn=\$(basename \$fna)
            prefix=\${bn%.fna.gz}
            gzip -cd \$fna > \${prefix}.fna
            prodigal -i \${prefix}.fna -f gff -o tmp.gff
            echo -e "##FASTA" | cat tmp.gff /dev/stdin \${prefix}.fna > gffs/\${prefix}.gff
        fi
    done
    panta.py -p init -g gffs/*.gff -o panta_refs -as -s -i 85 -c 20 -e 0.01 -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Pasa_Panta: 1.0
    END_VERSIONS
    """
}

process PASA {
    tag "$meta.id"
    label 'process_high'

    conda "python=3.7.10 biopython=1.79 bedtools=2.30.0 prodigal=2.6.3 cd-hit=4.8.1 blast=2.13.0 hmmer=3.3.2 diamond=2.0.14 mcl=14.137 mafft=7.526 parallel=20220222 numpy=1.21.6 scipy=1.7.3 numba=0.56.3 networkx=2.6.3 pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-adf032e06d3ecc4986e4c4efa160e481707acd5a:5d5dd40d2a8d9df574dfc603ccb059d9d3354e19-0' :
        'scgs/mulled-v2-adf032e06d3ecc4986e4c4efa160e481707acd5a:5d5dd40d2a8d9df574dfc603ccb059d9d3354e19-0' }"

    input:
    tuple val(meta), path(spades_out)
    path(panta_refs)

    output:
    tuple val(meta), path("${prefix}.scaffolds.fasta")          , emit: scaffolds
    tuple val(meta), path("${prefix}.ctg200.fasta")             , emit: ctg200
    tuple val(meta), path("${prefix}.ctgs.fasta")               , emit: ctg
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -arL ${panta_refs} panta_${prefix}
    prodigal -i ${spades_out}/contigs.fasta -f gff -o tmp.gff
    echo -e "##FASTA" | cat tmp.gff /dev/stdin ${spades_out}/contigs.fasta > ${prefix}.gff
    panta.py -p add -g ${prefix}.gff -o panta_${prefix} -as -s -i 85 -c 20 -e 0.01 -t ${task.cpus}
    pasa.py --data_dir panta_${prefix} --incomplete_sample_name ${prefix} --assem_dir $spades_out --output_fasta ${prefix}.scaffolds.fasta
    faFilterByLen.pl ${prefix}.scaffolds.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/_length.*\$//g' > ${prefix}.ctgs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Pasa_Panta: 1.0
    END_VERSIONS
    """
}
