nextflow.enable.types = true

process BBNORM {
    tag "$meta.id"
    label 'process_medium'

    conda "fastx_toolkit=0.0.14 bioconda::bbmap=39.01"
    container "scgs/mulled-v2-0f45a2e9949b9309cc37635f57bff7a66baf8095:86172d512030702a6bdb7b2cd7e301c3e1a14e56-1"

    input:
    tuple(meta: Map, reads: List<Path>)

    output:
    record(meta: meta, single_fastq: file('*_norm.fastq.gz', optional: true), fastq1: file('*_norm_R1.fastq.gz', optional: true), fastq2: file('*_norm_R2.fastq.gz', optional: true), log: file('*.log'), versions: file('versions.yml'))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = params.bulk ? "bulk" : "mda"
    def memory = ((task.memory.toGiga() * 0.8) as Integer) + 'g'
    if (meta.single_end) {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
        ln -s ${reads[0]} ${prefix}_norm.fastq.gz
    else
        bbnorm.sh in=${reads[0]} out=${prefix}_norm.fastq.gz $args threads=$task.cpus -Xmx${memory} &> ${prefix}.bbnorm.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
    } else {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
        ln -s ${reads[0]} ${prefix}_norm_R1.fastq.gz
        ln -s ${reads[1]} ${prefix}_norm_R2.fastq.gz
    else
        gzip -cd ${reads[0]} | fastx_renamer -n COUNT -i /dev/stdin -Q33 -z -o ${prefix}_rename_R1_fq.gz
        gzip -cd ${reads[1]} | fastx_renamer -n COUNT -i /dev/stdin -Q33 -z -o ${prefix}_rename_R2_fq.gz
        bbnorm.sh in=${prefix}_rename_R1_fq.gz in2=${prefix}_rename_R2_fq.gz \\
        out=${prefix}_norm_R1.fastq.gz out2=${prefix}_norm_R2.fastq.gz \\
        $args threads=$task.cpus -Xmx${memory} &> ${prefix}.bbnorm.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
    }
}
