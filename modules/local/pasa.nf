nextflow.enable.types = true

process PANTA {
    tag "panta"
    label 'process_high'

    conda "pypy=7.3.15 python=3.9.18 biopython=1.84 bedtools=2.30.0 prodigal=2.6.3 cd-hit=4.8.1 blast=2.13.0 hmmer=3.3.2 diamond=2.0.14 mcl=14.137 mafft=7.526 parallel=20220222 numpy=1.26.4 scipy=1.12.0 networkx=2.6.3 pandas=2.2.2 perl-bioperl=1.7.8"
    container "scgs/mulled-v2-073b771ca2dadccea705dbf1ddd01a7cf8acbd16:2dbb37a53c6b2b0022680b85f721d8d95f888d99-0"

    input:
    refs_fna: Bag<Path>

    output:
    record(db: file("panta_refs", type: "dir"))
    topic:
    file("versions.yml") >> 'local_versions'

    script:
    """
    mkdir -p gffs
    refs=(${refs_fna})
    for fna in \${refs[*]}; do
        if [[ \$fna == *.fna.gz ]]; then
            bn=\$(basename \$fna)
            prefix=\${bn%.fna.gz}
            gzip -cd \$fna > \${prefix}.fna
            prodigal -i \${prefix}.fna -c -m -q -f gff -o tmp.gff
            echo -e "##FASTA" | cat tmp.gff /dev/stdin \${prefix}.fna > gffs/\${prefix}.gff
        else
            bn=\$(basename \$fna)
            prefix=\${bn%.fna}
            prodigal -i \$fna -c -m -q -f gff -o tmp.gff
            echo -e "##FASTA" | cat tmp.gff /dev/stdin \$fna > gffs/\${prefix}.gff
        fi
    done
    panta.py -p init -g gffs/*.gff -o panta_refs -as -s -i 85 -c 50 -e 0.00001 -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Pasa_Panta: 1.0
    END_VERSIONS
    """
}

process PASA {
    tag "$meta.id"
    label 'process_high'

    conda "pypy=7.3.15 python=3.9.18 biopython=1.84 bedtools=2.30.0 prodigal=2.6.3 cd-hit=4.8.1 blast=2.13.0 hmmer=3.3.2 diamond=2.0.14 mcl=14.137 mafft=7.526 parallel=20220222 numpy=1.26.4 scipy=1.12.0 networkx=2.6.3 pandas=2.2.2 perl-bioperl=1.7.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-073b771ca2dadccea705dbf1ddd01a7cf8acbd16:2dbb37a53c6b2b0022680b85f721d8d95f888d99-0' :
        'scgs/mulled-v2-073b771ca2dadccea705dbf1ddd01a7cf8acbd16:2dbb37a53c6b2b0022680b85f721d8d95f888d99-0' }"

    input:
    tuple(meta: Map, spades_out: Path)
    panta_refs: Path

    output:
    record(meta: meta, scaffolds: file("*.scaffolds.fasta"), ctg200: file("*.pasa200.fasta"), ctg: file("*.pasa.fasta"))
    topic:
    file("versions.yml") >> 'local_versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -arL ${panta_refs} panta_${prefix}
    cp -arL $spades_out spades_for_pasa
    prodigal -i spades_for_pasa/contigs.fasta -c -m -q -f gff -o tmp.gff
    echo -e "##FASTA" | cat tmp.gff /dev/stdin spades_for_pasa/contigs.fasta > ${prefix}.gff
    panta.py -p add -g ${prefix}.gff -o panta_${prefix} -as -s -i 85 -c 50 -e 0.00001 -t ${task.cpus}
    pasa.py --data_dir panta_${prefix} --incomplete_sample_name ${prefix} --assem_dir spades_for_pasa --output_fasta ${prefix}.pasa.fasta
    fixSPAdesLen.py ${prefix}.pasa.fasta | sed 's/NODE_/PASA_/g' > ${prefix}.scaffolds.fasta
    faFilterByLen.pl ${prefix}.scaffolds.fasta 200 > ${prefix}.pasa200.fasta
    cat ${prefix}.pasa200.fasta | sed 's/_length.*\$//g' > ${prefix}.pasa.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Pasa_Panta: 1.0
    END_VERSIONS
    """
}
