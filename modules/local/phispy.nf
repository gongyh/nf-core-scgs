nextflow.enable.types = true

process PHISPY {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11.5 conda-forge::numpy=1.26.4 conda-forge::pandas=2.3.2 conda-forge::scikit-learn=1.7.1 bioconda::pyrodigal=3.6.3.post1 bioconda::phispy=4.2.21 conda-forge::pexpect=4.9.0"
    container "scgs/mulled-v2-429a3460971b0153ab4b5691b696eab3d551813d:54e9422a549b5e87e5486d5c5b9b5fcdfcca1bd7-0"

    input:
    tuple(meta: Map, gbk: Path)

    output:
    record(meta: meta, out_operon: file("${prefix}", type: 'dir'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    # prophages identification
    set +e
    PhiSpy.py ${gbk} -o ${prefix} --threads ${task.cpus} --color
    phispy_status=\$?
    set -e
    if [ "\$phispy_status" -eq 41 ]; then
        touch ${prefix}/NO_PROPHAGES_FOUND
    elif [ "\$phispy_status" -ne 0 ]; then
        exit "\$phispy_status"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phispy: 4.2.21
    END_VERSIONS
    """
}
