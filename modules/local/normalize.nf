process NORMALIZE {
    tag "${prefix}"

    input:
    path clean_reads
    val(single_end)

    output:
    path("*_norm*.fastq.gz"),                     emit: reads

    when:
    params.ass

    script:
    prefix = clean_reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?(\.bz2)?$/
    def R1 = clean_reads[0].toString()
    def mode = params.bulk ? "bulk" : "mda"
    if (single_end) {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
    ln -s $R1 ${prefix}_norm.fastq.gz
    else
    normalize-by-median.py -k 31 -C 40 --gzip -M 4e+9 -R ${prefix}_norm.report -o ${prefix}_norm.fastq.gz $R1
    fi
    """
    } else {
    def R2 = clean_reads[1].toString()
    """
    if [ \"${mode}\" == \"bulk\" ]; then
    ln -s $R1 ${prefix}_norm_R1.fastq.gz
    ln -s $R2 ${prefix}_norm_R2.fastq.gz
    else
    gzip -cd $R1 | fastx_renamer -n COUNT -i /dev/stdin -Q33 -z -o ${prefix}_rename_R1_fq.gz
    gzip -cd $R2 | fastx_renamer -n COUNT -i /dev/stdin -Q33 -z -o ${prefix}_rename_R2_fq.gz
    interleave-reads.py ${prefix}_rename_R1_fq.gz ${prefix}_rename_R2_fq.gz | normalize-by-median.py -k 31 -C 40 -M 4e+9 -p --gzip -R ${prefix}_norm.report -o ${prefix}_nbm.fastq.gz /dev/stdin
    split-paired-reads.py -1 ${prefix}_norm_R1.fastq.gz -2 ${prefix}_norm_R2.fastq.gz --gzip ${prefix}_nbm.fastq.gz
    fi
    """
    }
}
