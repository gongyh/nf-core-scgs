process PANTA_INIT {
    tag "panta_refs"
    label 'process_high'

    conda "python=3.7.10 biopython=1.79 bedtools=2.30.0 prodigal=2.6.3 cd-hit=4.8.1 blast=2.13.0 hmmer=3.3.2 diamond=2.0.14 mcl=14.137 mafft=7.526 parallel=20220222 numpy=1.21.6 scipy=1.7.3 numba=0.56.3 networkx=2.6.3 pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-adf032e06d3ecc4986e4c4efa160e481707acd5a:5d5dd40d2a8d9df574dfc603ccb059d9d3354e19-0' :
        'scgs/mulled-v2-adf032e06d3ecc4986e4c4efa160e481707acd5a:5d5dd40d2a8d9df574dfc603ccb059d9d3354e19-0' }"

    input:
    path(refs)

    output:
    path("panta_refs")                , emit: out
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    panta.py -p init -a ${refs} -o panta_refs -as -s -i 85 -c 20 -e 0.01 -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Pasa_Panta: 1.0
    END_VERSIONS
    """
}

process PANTA_ADD {
    tag "$meta.id"
    label 'process_high'

    conda "python=3.7.10 biopython=1.79 bedtools=2.30.0 prodigal=2.6.3 cd-hit=4.8.1 blast=2.13.0 hmmer=3.3.2 diamond=2.0.14 mcl=14.137 mafft=7.526 parallel=20220222 numpy=1.21.6 scipy=1.7.3 numba=0.56.3 networkx=2.6.3 pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-adf032e06d3ecc4986e4c4efa160e481707acd5a:5d5dd40d2a8d9df574dfc603ccb059d9d3354e19-0' :
        'scgs/mulled-v2-adf032e06d3ecc4986e4c4efa160e481707acd5a:5d5dd40d2a8d9df574dfc603ccb059d9d3354e19-0' }"

    input:
    tuple val(meta), path(contigs)
    path(panta_refs)

    output:
    path("panta_${prefix}")             , emit: out
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -ar ${panta_refs} panta_${prefix}
    panta.py -p add -a $contigs -o panta_${prefix} -as -s -i 85 -c 20 -e 0.01 -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Pasa_Panta: 1.0
    END_VERSIONS
    """
}
