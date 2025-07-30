process QUICKMERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "scgs::quickmerge-vt=0.4 bioconda::seqkit=2.10.0 conda-forge::biopython=1.85 bioconda::perl-bioperl=1.7.8 bioconda::seqtk=1.4"
    container "scgs/mulled-v2-3d7dbca3694e0bc412a351a0d57cbefcab830270:13f4c4f6b30adc2e47d2dfca9753708b5d6e5cc7-0"

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
        merge_wrapper.py -v -t ${task.cpus} -l 0 -pre ${prefix} -ml 100 \${refass_contig} tmp.fasta
        cp -f merged_${prefix}.fasta tmp.fasta
    done
    # clean up read id
    seqtk rename merged_${prefix}.fasta ${prefix}_ | sed 's/ .*\$//g' > ${prefix}.hybrid.fasta
    # remove short contigs
    faFilterByLen.pl ${prefix}.hybrid.fasta 200 > ${prefix}.hybrid200.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quickmerge: 0.4
    END_VERSIONS
    """
}
