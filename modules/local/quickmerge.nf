process QUICKMERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::quickmerge=0.3 bioconda::seqkit=2.10.0 conda-forge::biopython=1.85 bioconda::perl-bioperl=1.7.8"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-446fe851503c98aa856d04200e6a396c80d50643:181987cf2fc30068e34f3bcc5464c8219c8d6ec2-1'
        : 'scgs/mulled-v2-446fe851503c98aa856d04200e6a396c80d50643:181987cf2fc30068e34f3bcc5464c8219c8d6ec2-1'}"

    input:
    tuple val(meta), path(denovo_contigs) // denovo assembled assembly
    tuple val(meta), path(refass_contigs) // reference guided assembly, after scaffolding

    output:
    tuple val(meta), path("${prefix}_merged200.fasta"),   emit: merged_assembly
    tuple val(meta), path("${prefix}_clean.fasta"),       emit: merged_clean
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ## merge denovo and ref-based SAGs
    merge_wrapper.py -pre ${prefix} -ml 100 ${denovo_contigs} ${refass_contigs}
    # append unaligned seqs
    cut -f1 aln_summary_${prefix}.tsv | grep -v REF | sort | uniq > aln_${prefix}.ids
    seqkit grep -v -n -f aln_${prefix}.ids ${refass_contigs} > unaln_${prefix}.fasta
    cat merged_${prefix}.fasta unaln_${prefix}.fasta > merged2_${prefix}.fasta
    # extract contigs from scaffolds
    scf2ctg.py merged2_${prefix}.fasta ${prefix}_merged.fasta
    # remove short contigs
    faFilterByLen.pl ${prefix}_merged.fasta 200 > ${prefix}_merged200.fasta
    cat ${prefix}_merged200.fasta | sed 's/_length.*\$//g' > ${prefix}_clean.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quickmerge: 0.3
    END_VERSIONS
    """
}
