process NORMALIZE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_norm*.fastq.gz"),  emit: reads

    when:
    params.ass

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def mode = params.bulk ? "bulk" : "mda"
    if (meta.single_end) {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
    ln -s ${reads[0]} ${prefix}_norm.fastq.gz
    else
    normalize-by-median.py -k 31 -C 40 --gzip -M 4e+9 -R ${prefix}_norm.report -o ${prefix}_norm.fastq.gz $R1
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
    interleave-reads.py ${prefix}_rename_R1_fq.gz ${prefix}_rename_R2_fq.gz | normalize-by-median.py -k 31 -C 40 -M 4e+9 -p --gzip -R ${prefix}_norm.report -o ${prefix}_nbm.fastq.gz /dev/stdin
    split-paired-reads.py -1 ${prefix}_norm_R1.fastq.gz -2 ${prefix}_norm_R2.fastq.gz --gzip ${prefix}_nbm.fastq.gz
    fi
    """
    }
}