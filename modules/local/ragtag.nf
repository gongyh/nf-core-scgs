process RAGTAG {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::ragtag=2.1.0 bioconda::seqkit=2.10.0"
    container "scgs/mulled-v2-561a2673ebd796b3ddd2822d3f38440d215223c5:6b65b2e5d7cc53084c2dc5fec2260d8adbee49f1-0"

    input:
    tuple val(meta), path(refass_contigs), path(denovo_contigs) // ref and denovo assemblies
    path(refs_fna)

    output:
    tuple val(meta), path("${prefix}_scaffolds.fasta"),    emit: scaffolded_assembly
    tuple val(meta), path("${prefix}.denovo.clean.fasta"), emit: denovo_assembly
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
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
        ragtag.py scaffold -f 200 -o ref_\${pre} -u -t ${task.cpus} --mm2-params '-x asm10' \${pre}.fna ${refass_contigs}
        minimap2 -x asm20 -t ${task.cpus} -o \${pre}_mapped.paf \${pre}.fna ${denovo_contigs}
    done
    ragtag.py merge -l 200 -u -o ragtag_merge ${refass_contigs} ref_*/*.agp
    cat *_mapped.paf | cut -f1 | sort | uniq > denovo_clean.gids
    seqkit grep -f denovo_clean.gids -o ${prefix}.denovo.clean.fasta ${denovo_contigs}
    ragtag.py patch -f 200 -o ragtag_patch -t ${task.cpus} --nucmer-params '--maxmatch -l 50 -c 250' -u --fill-only ragtag_merge/ragtag.merge.fasta ${prefix}.denovo.clean.fasta
    ragtag.py correct -f 200 -o ragtag_correct -u -t ${task.cpus} --mm2-params '-x asm10' --intra ${prefix}.denovo.clean.fasta ragtag_patch/ragtag.patch.fasta

    cp ragtag_correct/ragtag.correct.fasta ${prefix}_scaffolds.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RagTag: \$(echo \$(ragtag.py -v | sed 's/v//'))
    END_VERSIONS
    """
}
