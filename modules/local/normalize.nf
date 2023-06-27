process NORMALIZE {
    tag "${meta.id}"

    conda "bioconda::khmer=3.0.0a3 bioconda::fastx_toolkit=0.0.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-afeccb6637ecc3e429a8f7f6e6713be70eff3d40:ee0ebbe2f959481c603cf90cda1c2026613505ef-0' :
        'scgs/mulled-v2-afeccb6637ecc3e429a8f7f6e6713be70eff3d40:ee0ebbe2f959481c603cf90cda1c2026613505ef-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_norm*.fastq.gz"),  emit: reads

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = params.bulk ? "bulk" : "mda"
    if (meta.single_end) {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
        ln -s ${reads[0]} ${prefix}_norm.fastq.gz
    else
        normalize-by-median.py -k 31 -C 40 --gzip -M ${task.memory.toGiga()}G -R ${prefix}_norm.report -o ${prefix}_norm.fastq.gz $R1
    fi
    """
    } else {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
        ln -s ${reads[0]} ${prefix}_norm_R1.fastq.gz
        ln -s ${reads[1]} ${prefix}_norm_R2.fastq.gz
    else
        gzip -cd ${reads[0]} | fastx_renamer -n COUNT -i /dev/stdin -Q33 -z -o ${prefix}_rename_R1_fq.gz
        gzip -cd ${reads[1]} | fastx_renamer -n COUNT -i /dev/stdin -Q33 -z -o ${prefix}_rename_R2_fq.gz
        interleave-reads.py ${prefix}_rename_R1_fq.gz ${prefix}_rename_R2_fq.gz | normalize-by-median.py -k 31 -C 40 -M ${task.memory.toGiga()}G -p --gzip -R ${prefix}_norm.report -o ${prefix}_nbm.fastq.gz /dev/stdin
        split-paired-reads.py -1 ${prefix}_norm_R1.fastq.gz -2 ${prefix}_norm_R2.fastq.gz --gzip ${prefix}_nbm.fastq.gz
    fi
    """
    }
}
