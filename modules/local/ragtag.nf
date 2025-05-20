process RAGTAG {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::ragtag=2.1.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ragtag:2.1.0--pyhb7b1952_0'
        : 'biocontainers/ragtag:2.1.0--pyhb7b1952_0'}"

    input:
    tuple val(meta), path(refass_contigs) // reference guided assembly
    tuple val(meta2), path(denovo_contigs) // denovo assembled assembly
    tuple path(refs_fna)

    output:
    tuple val(meta), path("${prefix}_scaffolds.fasta"),   emit: scaffolded_assembly
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    refs=(${refs_fna})
    for fna in \${refs[*]}; do
        bn=\$(basename \$fna)
        if [[ \$fna == *.fna.gz ]]; then
            pre=\${bn%.fna.gz}
            gzip -cd \$fna > \${pre}.fna
        else
            pre=\${bn%.fna}
        fi
        ragtag.py scaffold -f 200 -o ref_\${pre} -u -t ${task.cpus} \${pre}.fna ${refass_contigs}
    done
    ragtag.py merge -l 200 -u -o ragtag_merge ${refass_contigs} ref_*/*.agp
    ragtag.py patch -f 200 -o ragtag_patch -t ${task.cpus} -u --fill-only ragtag_merge/ragtag.merge.fasta ${denovo_contigs}
    ragtag.py correct -f 200 -o ragtag_correct -u -t ${task.cpus} --intra ${denovo_contigs} ragtag_patch/ragtag.patch.fasta
    
    cp ragtag_correct/ragtag.correct.fasta ${prefix}_scaffolds.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ragtag: \$(echo \$(ragtag.py -v | sed 's/v//'))
    END_VERSIONS
    """
}
