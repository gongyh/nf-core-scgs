nextflow.enable.types = true

process PROMPREDICT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.11.5 conda-forge::numpy=1.26.4 conda-forge::pandas=2.3.2 conda-forge::scikit-learn=1.7.1 bioconda::pyrodigal=3.6.3.post1 bioconda::phispy=4.2.21 conda-forge::pexpect=4.9.0"
    container "scgs/mulled-v2-429a3460971b0153ab4b5691b696eab3d551813d:54e9422a549b5e87e5486d5c5b9b5fcdfcca1bd7-0"

    input:
    tuple(meta: Map, genome_fasta: Path)

    output:
    record(meta: meta, out_operon: file("${prefix}", type: 'dir'), versions: file('versions.yml'))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp ${genome_fasta} ${prefix}.fna
    mkdir -p ${prefix}
    # promoter identification
    PromPredict_genome_V1.py --genome_fasta ${prefix}.fna
    cp ${prefix}_*.txt ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PromPredict: V1
    END_VERSIONS
    """
}
