process PREPARE_REASSEMBLY_TARGETS {
    tag "prepare_targets"
    label 'process_low'
    container "community.wave.seqera.io/library/dnaberts:7a7299083f265248"

    input:
    path checkm2_report
    path taxvamb_report
    path bins_dir
    path script_py
    output:
    path "targets", emit: targets_dir

    script:
    """
    python3 ${script_py} \
        --checkm2 ${checkm2_report} \
        --taxonomy ${taxvamb_report} \
        --bindir ${bins_dir} \
        --outdir targets
    """
}

process POST_BINNING_REASSEMBLY {
    tag "reassemble"
    label 'process_high'
    container "community.wave.seqera.io/library/bowtie2_samtools_spades:cdb9830318b544b0" 

    input:
    path targets_dir
    path all_reads

    output:
    path "reassembled_bins", emit: reassembled

    script:
    """
    mkdir -p reassembled_bins
    if ! ls ${targets_dir}/*.fa >/dev/null 2>&1; then
        echo "错误：没有找到任何目标 .fa 文件！"
        exit 1
    fi

    for target in ${targets_dir}/*.fa; do
        name=\$(basename \$target .fa)
        echo "======================================"
        echo "开始重组装靶标: \$name"
        cp -L \$target ./target.fa
        bowtie2-build ./target.fa ./\${name}_db > /dev/null 2>&1
        > \${name}_mapped_1.fq
        > \${name}_mapped_2.fq
        for r1 in ${all_reads}; do
            if [[ "\$r1" == *_1.fastq ]]; then
                r2=\${r1/_1.fastq/_2.fastq}

                bowtie2 -x ./\${name}_db -1 \$r1 -2 \$r2 -S tmp.sam -p 4
                samtools view -b -F 12 tmp.sam | samtools fastq -1 tmp_1.fq -2 tmp_2.fq -0 /dev/null -s /dev/null -
                cat tmp_1.fq >> \${name}_mapped_1.fq
                cat tmp_2.fq >> \${name}_mapped_2.fq
                rm -f tmp.sam tmp_1.fq tmp_2.fq
            fi
        done
        spades.py --careful -1 ./\${name}_mapped_1.fq -2 ./\${name}_mapped_2.fq -o reassembled_bins/\${name}_out --threads 8
        rm -rf ./target.fa ./\${name}_db ./\${name}_mapped_1.fq ./\${name}_mapped_2.fq
    done
    """
}
