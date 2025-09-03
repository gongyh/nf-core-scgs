process PHISPY {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11.5 conda-forge::numpy=1.26.4 conda-forge::pandas=2.3.2 conda-forge::scikit-learn=1.7.1 bioconda::pyrodigal=3.6.3.post1 bioconda::phispy=4.2.21 conda-forge::pexpect=4.9.0"
    container "scgs/mulled-v2-429a3460971b0153ab4b5691b696eab3d551813d:54e9422a549b5e87e5486d5c5b9b5fcdfcca1bd7-0"

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("$prefix"), emit: out_operon
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    # prophages identification
    PhiSpy.py ${gbk} -o ${prefix} --threads ${task.cpus} --color

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phispy: 4.2.21
    END_VERSIONS
    """
}
