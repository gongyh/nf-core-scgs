process SPADES_JOINT {
    tag "joint"
    label 'process_high'          // 需要大内存（如 64GB）
    publishDir "${params.outdir}/joint", mode: 'copy'

    input:
    path r1
    path r2
    path s

    output:
    path "joint_assembly/contigs.fasta", emit: contigs
    path "joint_assembly/", emit: assembly_dir
    path "versions.yml", emit: versions

    script:
    """
    spades.py --sc --careful -k 55,77,99 -t ${task.cpus} -m ${task.memory.toGiga()} \
        -1 ${r1} -2 ${r2} -s ${s} -o joint_assembly

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//')
    END_VERSIONS
    """
}
