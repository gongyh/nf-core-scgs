process QUICKMERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::quickmerge=0.3 bioconda::seqkit=2.10.0 conda-forge::biopython=1.85 bioconda::perl-bioperl=1.7.8 bioconda::seqtk=1.4"
    container "scgs/mulled-v2-d417af7602b66a7a02bee82c7dd6399da6f61ce0:d831d87d4fdb108118b1d07ed3b32621cd2472f2-0"

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
    for refass_contig in \${refass_contigs[*]}; do
        merge_wrapper.py -pre ${prefix} -ml 100 \${refass_contig} tmp.fasta
        cp -f merged_${prefix}.fasta tmp.fasta
    done
    # clean up read id
    seqtk rename merged_${prefix}.fasta ${prefix}_ | sed 's/ .*\$//g' > ${prefix}.hybrid.fasta
    # remove short contigs
    faFilterByLen.pl ${prefix}.hybrid.fasta 200 > ${prefix}.hybrid200.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quickmerge: 0.3
    END_VERSIONS
    """
}
