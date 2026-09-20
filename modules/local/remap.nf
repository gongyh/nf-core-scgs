nextflow.enable.types = true

process REMAP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bowtie2=2.4.4 bioconda::samtools=1.17"
    container "scgs/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:7a723346025146e191fdbf519b8131c258b9eeab-0"

    input:
    tuple(meta: Map, reads: List<Path>, index: Path)
    allow_multi_align: Boolean

    output:
    record(meta: meta, bam: file("*_ass.sort.bam"), bai: file("*_ass.sort.bam.bai"), mqc_tsv: file("remap_mqc.tsv"))
    topic:
    file("versions.yml") >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def index_name = meta.id_index ? "${meta.id_index}" : prefix
    def filtering = allow_multi_align ? '' : "| samtools view -b -q 40 -F 4 -F 256 -"
    if (meta.single_end) {
    """
    bowtie2 -x ${index_name}Bowtie2Index/${index_name} -p ${task.cpus} -U ${reads} 2> bowtie2.log | samtools view -bT ${index_name}Bowtie2Index - $filtering > ${prefix}_ass.bam
    samtools sort -o ${prefix}_ass.sort.bam ${prefix}_ass.bam
    samtools index ${prefix}_ass.sort.bam

    ALIGN_RATE=\$(grep "overall alignment rate" bowtie2.log | awk '{print \$1}' | tr -d '%')
    printf "Metric\tValue\n" > remap_mqc.tsv
    printf "Overall alignment rate (%%)\t\${ALIGN_RATE}\n" >> remap_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
    """
    } else {
    """
    bowtie2 --no-mixed --no-discordant -X 1000 -x ${index_name}Bowtie2Index/${index_name} -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} 2> bowtie2.log | samtools view -bT ${index_name}Bowtie2Index - $filtering > ${prefix}_ass.bam
    samtools sort -o ${prefix}_ass.sort.bam ${prefix}_ass.bam
    samtools index ${prefix}_ass.sort.bam

    ALIGN_RATE=\$(grep "overall alignment rate" bowtie2.log | awk '{print \$1}' | tr -d '%')
    printf "Metric\tValue\n" > remap_mqc.tsv
    printf "Overall alignment rate (%%)\t\${ALIGN_RATE}\n" >> remap_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
    """
    }
}
