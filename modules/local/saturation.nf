process SATURATION {
    tag { $prefix }
    publishDir path: "${params.outdir}/saturation", mode: 'copy'

    input:
    path(reads)
    val(single_end)

    output:
    path("${prefix}_kmer.pdf")
    path("${prefix}_cov31_*.csv")

    when:
    params.saturation

    script:
    prefix = reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    R1 = reads[0].toString()
    if (single_end) {
    """
    fastp -i $R1 -A -G -Q -L -s 10 -d 0 -o ${prefix}_split.fq.gz
    for i in {1..10}; do
    mccortex31 build --kmer 31 --sample \$i -t ${task.cpus} -Q 20 -m 8G \
        --seq \${i}.${prefix}_split.fq.gz \${i}.k31.ctx
    if [ \$i == 1 ]; then
        mccortex31 clean -t ${task.cpus} -m 8G -U10 -T16 -f -o null -C ${prefix}_cov31_p\${i}.csv 0:\${i}.k31.ctx
        cp -f 1.k31.ctx tmp_clean31.ctx
    else
        mccortex31 join -m 8G --out merged_clean31.ctx 0:\${i}.k31.ctx 0:tmp_clean31.ctx
        mccortex31 clean -t ${task.cpus} -m 8G -U10 -T16 -f -o null -C ${prefix}_cov31_p\${i}.csv 0:merged_clean31.ctx
        mv -f merged_clean31.ctx tmp_clean31.ctx
    fi
    done
    KmerDensity.R \$PWD ${prefix}
    """
    } else {
    R2 = reads[1].toString()
    """
    fastp -i $R1 -I $R2 -A -G -Q -L -s 10 -d 0 -o ${prefix}_split_R1.fq.gz -O ${prefix}_split_R2.fq.gz
    for i in {1..10}; do
    mccortex31 build --kmer 31 --sample \$i -t ${task.cpus} -Q 20 -m 8G \
        --seq2 \${i}.${prefix}_split_R1.fq.gz:\${i}.${prefix}_split_R2.fq.gz \${i}.k31.ctx
    if [ \$i == 1 ]; then
        mccortex31 clean -t ${task.cpus} -m 8G -U10 -T16 -f -o null -C ${prefix}_cov31_p\${i}.csv 0:\${i}.k31.ctx
        cp -f 1.k31.ctx tmp_clean31.ctx
    else
        mccortex31 join -m 8G --out merged_clean31.ctx 0:\${i}.k31.ctx 0:tmp_clean31.ctx
        mccortex31 clean -t ${task.cpus} -m 8G -U10 -T16 -f -o null -C ${prefix}_cov31_p\${i}.csv 0:merged_clean31.ctx
        mv -f merged_clean31.ctx tmp_clean31.ctx
    fi
    done
    KmerDensity.R \$PWD ${prefix}
    """
    }

}
