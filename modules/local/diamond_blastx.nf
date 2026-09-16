nextflow.enable.types = true

process DIAMOND_BLASTX {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::diamond=2.0.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_0' :
        'biocontainers/diamond:2.0.15--hb97b32f_0' }"

    input:
    tuple(meta: Map, contigs: Path)
    tuple(_nt_meta: Object, nt_out: Path)
    uniprot: Path
    uniprot_taxids: Path
    has_uniprot: Boolean

    output:
    record(meta: meta, uniprot: file("*_uniprot.taxified.out"), contigs: file("*.fasta"), nt: file("*.out"), has_uniprot: has_uniprot, versions: file("versions.yml"), out_put: file("*_uniprot.*"))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!has_uniprot) {
    """
    touch ${prefix}_uniprot.out
    touch ${prefix}_uniprot.taxified.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(echo \$(diamond version 2>&1) | sed 's/^.*diamond version //; s/Using.*\$//')
    END_VERSIONS
    """
    } else {
    """
    diamond blastx --query $contigs --db $uniprot -p ${task.cpus} -o ${prefix}_uniprot.out \
        --outfmt 6 --sensitive --max-target-seqs 1 --evalue ${params.evalue} -b ${params.blockSize}
    blobtools taxify -f ${prefix}_uniprot.out -m ${uniprot_taxids} -s 0 -t 2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(echo \$(diamond version 2>&1) | sed 's/^.*diamond version //; s/Using.*\$//')
    END_VERSIONS
    """
    }
}
