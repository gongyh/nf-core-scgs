process QUICKMERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::quickmerge=0.3 bioconda::seqkit=2.10.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-71a73be5aed9650e7416aa6b5810b891b17cfdce:c342bdf56a62c7f1a18890e7a4504765e808b6c8-0'
        : 'scgs/mulled-v2-71a73be5aed9650e7416aa6b5810b891b17cfdce:c342bdf56a62c7f1a18890e7a4504765e808b6c8-0'}"

    input:
    tuple val(meta), path(denovo_contigs) // denovo assembled assembly
    tuple val(meta), path(refass_contigs) // reference guided assembly, after scaffolding

    output:
    tuple val(meta), path("${prefix}_merged.fasta"),   emit: merged_assembly
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ## merge denovo and ref-based SAGs
    merge_wrapper.py -pre ${prefix} -ml 100 ${denovo_contigs} ${refass_contigs}
    # append unaligned seqs
    cut -f1 aln_summary_${prefix}.tsv | grep -v REF | sort | uniq > aln_${prefix}.ids
    seqkit grep -v -n -f aln_${prefix}.ids ${refass_contigs} > unaln_${prefix}.fasta
    cat merged_${prefix}.fasta unaln_${prefix}.fasta > merged2_${prefix}.fasta
    # extract contigs from scaffolds
    python scf2ctg.py merged2_${prefix}.fasta ${prefix}_merged.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quickmerge: 0.3
    END_VERSIONS
    """
}
