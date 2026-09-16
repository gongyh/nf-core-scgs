nextflow.enable.types = true

process PROKKA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::prokka=1.14.6 bioconda::bedops=2.4.38"
    container "scgs/mulled-v2-1e40df84b5b2d0a934c357a759500c269d2eb793:81460e1910925aa1427c823417f44d2739507564-0"

    input:
    tuple(meta: Map, contigs: Path)
    proteins: List<Path>

    output:
    record(meta: meta, prokka_for_split: file("*", type: "dir"), faa: file("*/*.faa"), gbk: file("*/*.gbk"), versions: file("versions.yml"))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    """
    prokka --outdir $prefix --prefix $prefix --strain $prefix --addgenes --addmrna --cpus ${task.cpus} $proteins_opt $contigs
    sed '/^##FASTA/Q' ${prefix}/${prefix}.gff > ${prefix}/${prefix}_noseq.gff
    gff2bed < ${prefix}/${prefix}_noseq.gff | cut -f1,4 | grep _gene | sed 's/_gene//g' > ${prefix}/${prefix}_ctg_genes.tsv
    prokka_postprocess.py ${prefix}/${prefix}_ctg_genes.tsv ${prefix}/${prefix}.tsv > ${prefix}/${prefix}_all.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka -v 2>&1) | sed 's/^.*prokka //; s/Using.*\$//')
    END_VERSIONS
    """
}
