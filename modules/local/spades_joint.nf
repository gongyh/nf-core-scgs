process SPADES {
    tag "joint"
    label 'process_high'
    publishDir "${params.outdir}/joint", mode: 'copy'

    input:
    path p1_files   
    path p2_files   
    path s_files    

    output:
    path "super_contigs.contigs.fa.gz", emit: contigs
    path "versions.yml", emit: versions

    script:
    """
    cat ${p1_files.join(' ')} > all_R1.fastq.gz
    cat ${p2_files.join(' ')} > all_R2.fastq.gz
    cat ${s_files.join(' ')}  > all_S.fastq.gz

    spades.py --sc --careful -k 55,77,99 -t ${task.cpus} -m ${task.memory.toGiga()} \\
        -1 all_R1.fastq.gz -2 all_R2.fastq.gz -s all_S.fastq.gz -o joint_assembly

    gzip -c joint_assembly/contigs.fasta > super_contigs.contigs.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/Using.*\$//')
    END_VERSIONS
    """
}
