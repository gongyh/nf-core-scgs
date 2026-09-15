process QUICKMERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "scgs::quickmerge-vt=0.4 bioconda::ragtag=2.1.0 bioconda::seqkit=2.10.0 conda-forge::biopython=1.85 bioconda::perl-bioperl=1.7.8 bioconda::seqtk=1.4"
    container "scgs/mulled-v2-3c99dbe67a0d01cc10a223e8f82778c618460187:2d2622edca5a7d6580a2ba583efd2d06d593b784-0"

    input:
    tuple val(meta), path(denovo_contigs), path(refass_contigs) // denovo and ref-guided assembled assemblies

    output:
    tuple val(meta), path("${prefix}.hybrid200.fasta"),   emit: merged_assembly
    tuple val(meta), path("${prefix}.hybrid.fasta"),      emit: merged_clean
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ## merge denovo and ref-based SAGs
    cp ${denovo_contigs} tmp.fasta
    refass_contigs=(${refass_contigs})
    # merge using QuickMerge
    for refass_contig in \${refass_contigs[*]}; do
        # correct ref-based assembly using RagTag
        ragtag.py correct -f 100 -b 50 -o ragtag_correct -u -t ${task.cpus} --mm2-params '-x asm10' --intra ${denovo_contigs} \${refass_contig}
        cp ragtag_correct/ragtag.correct.fasta corrected.\${refass_contig} # break intra missasseblies
        rm -rf ragtag_correct
        # merge
        merge_wrapper.py -v -t ${task.cpus} -l 1000 -pre ${prefix} -ml 200 tmp.fasta corrected.\${refass_contig}
        cp -f merged_${prefix}.fasta tmp.fasta
    done
    # clean up read id
    seqtk rename merged_${prefix}.fasta ${prefix}_ | sed 's/ .*\$//g' > ${prefix}.hybrid.fasta
    # remove short contigs
    faFilterByLen.pl ${prefix}.hybrid.fasta 200 > ${prefix}.hybrid200.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quickmerge: 0.4
        RagTag: \$(echo \$(ragtag.py -v | sed 's/v//'))
    END_VERSIONS
    """
}
