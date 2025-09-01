process UNIOP {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11.5 conda-forge::numpy=1.26.4 conda-forge::pandas=2.3.2 conda-forge::scikit-learn=1.7.1 bioconda::prodigal=2.6.3"
    container "scgs/mulled-v2-55a8e60d9443700fc30f7dcbd9e818915628799a:18ac8636f36b21579062b5b7aebc003241570d3f-0"

    input:
    tuple val(meta), path(prodigal_faa)

    output:
    tuple val(meta), path("$prefix"), emit: out_operon
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    UniOP.py -a ${prodigal_faa} -t ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UniOP: 1.0
    END_VERSIONS
    """
}
