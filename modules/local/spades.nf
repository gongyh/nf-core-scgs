nextflow.enable.types = true

process SPADES {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::spades=3.15.5 bioconda::perl-bioperl=1.7.8 conda-forge::python=3.10.14"
    container "scgs/mulled-v2-5524a20c8f39de906b127a66052c67b51c9a9ce1:c8e22953d04dee6a4da05f7a131bbd081ad78651-0"

    input:
    tuple(meta: Map, reads: List<Path>)

    output:
    record(meta: meta, corrected_read: file("*.corrected_R1.fastq.gz", optional: true), corrected_read2: file("*.corrected_R2.fastq.gz", optional: true), contig: file("*.contigs.fasta"), contig_path: file("*.contigs.paths"), contig_graph: file("*.spades_out/*.contigs.gfa"), ctg200: file("*.ctg200.fasta"), ctg: file("*.ctgs.fasta"), assembly: file("*.spades_out", type: "dir"), mqc_tsv: file("spades_joint_mqc.tsv"))
    topic:
    file("versions.yml") >> 'local_versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def mode = params.bulk ? "--cov-cutoff auto --careful" : "--sc --careful"
    mode = params.mg ? "--meta" : "--sc --careful"
    def rcl = meta.single_end ? "-s ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def se_pe = meta.single_end ? "SE" : "PE"
    """
    spades.py ${rcl} ${mode} ${args} -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out

    corrected_files=( \$(ls ${prefix}.spades_out/corrected/*.cor.fastq.gz 2>/dev/null | grep -v unpaired) )
    cp "\${corrected_files[0]}" ${prefix}.corrected_R1.fastq.gz
    if [ "${se_pe}" = "PE" ]; then
        cp "\${corrected_files[1]}" ${prefix}.corrected_R2.fastq.gz
    fi

    if [ "${args}" = "--only-error-correction" ]; then
        touch ${prefix}.contigs.fasta ${prefix}.contigs.paths ${prefix}.spades_out/${prefix}.contigs.gfa
        touch ${prefix}.ctg200.fasta ${prefix}.ctgs.fasta
    else
        cp ${prefix}.spades_out/assembly_graph_after_simplification.gfa ${prefix}.spades_out/${prefix}.contigs.gfa
        cp ${prefix}.spades_out/contigs.paths ${prefix}.spades_out/contigs.paths_raw
        cp ${prefix}.spades_out/scaffolds.paths ${prefix}.spades_out/scaffolds.paths_raw
        correctPaths.py ${prefix}.spades_out/contigs.fasta ${prefix}.spades_out/contigs.paths_raw ${prefix}.spades_out/contigs.paths
        correctPaths.py ${prefix}.spades_out/scaffolds.fasta ${prefix}.spades_out/scaffolds.paths_raw ${prefix}.spades_out/scaffolds.paths
        cp ${prefix}.spades_out/contigs.paths ${prefix}.contigs.paths
        cp ${prefix}.spades_out/contigs.fasta ${prefix}.contigs.fasta
        faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
        cat ${prefix}.ctg200.fasta | sed 's/_length.*\$//g' > ${prefix}.ctgs.fasta
    fi

    CONTIGS="${prefix}.contigs.fasta"
    if [ -f "\$CONTIGS" ]; then
        awk '/^>/ {if (seqlen) print seqlen; seqlen=0; next} { seqlen += length(\$0) } END {if (seqlen) print seqlen}' "\$CONTIGS" | sort -rn > lengths.txt
        TOTAL=\$(awk '{sum+=\$1} END {print sum}' lengths.txt)
        NUM=\$(wc -l < lengths.txt)
        LONGEST=\$(head -1 lengths.txt)
        half=\$((TOTAL / 2))
        cum=0; N50=0
        while read len; do
            cum=\$((cum + len))
            if [ \$cum -ge \$half ]; then
                N50=\$len
                break
            fi
        done < lengths.txt
        rm lengths.txt
    else
        NUM=0; TOTAL=0; N50=0; LONGEST=0
    fi

    printf "Metric\tValue\n" > spades_joint_mqc.tsv
    printf "Number of contigs (>=200bp)\t\${NUM}\n" >> spades_joint_mqc.tsv
    printf "Total assembly size (bp)\t\${TOTAL}\n" >> spades_joint_mqc.tsv
    printf "N50 (bp)\t\${N50}\n" >> spades_joint_mqc.tsv
    printf "Longest contig (bp)\t\${LONGEST}\n" >> spades_joint_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    spades: \$(echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//; s/Using.*\$//')
END_VERSIONS
    """
}
