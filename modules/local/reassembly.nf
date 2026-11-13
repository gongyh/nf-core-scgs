nextflow.enable.types = true

process PREPARE_REASSEMBLY_TARGETS {
    tag "prepare_targets"
    label 'process_low'
    container "community.wave.seqera.io/library/dnaberts:7a7299083f265248"

    input:
    checkm2_report: Path
    taxvamb_report: Path
    bins_dir: Path
    output:
    targets_dir = file("targets", type: 'dir')

    script:
    def script_path = "${projectDir}/bin/prepare_reassembly_targets.py"
    """
    python3 ${script_path} \\
        --checkm2 ${checkm2_report} \\
        --taxonomy ${taxvamb_report} \\
        --bindir ${bins_dir} \\
        --outdir targets
    """
}

process POST_BINNING_REASSEMBLY {
    tag "reassemble"
    label 'process_high'
    container "community.wave.seqera.io/library/bowtie2_samtools_spades:cdb9830318b544b0"

    input:
    targets_dir: Path
    all_reads: List<Path>

    output:
    reassembled_bins = file("reassembled_bins", type: 'dir')

    script:
    """
    mkdir -p reassembled_bins
    if ! ls ${targets_dir}/*.fa >/dev/null 2>&1; then
        echo "Error: No target .fa files found!"
        exit 1
    fi

    for target in ${targets_dir}/*.fa; do
        name=\$(basename \$target .fa)
        echo "======================================"
        echo "Starting reassembly for target: \$name"
        cp -L \$target ./target.fa
        bowtie2-build ./target.fa ./\${name}_db > /dev/null 2>&1
        > \${name}_mapped_1.fq
        > \${name}_mapped_2.fq
        > \${name}_mapped_single.fq
        for fq in ${all_reads.join(' ')}; do
            if [[ "\$fq" == *_R1.fastq* || "\$fq" == *_1.fastq* ]]; then
                r1=\$fq
                r2=\${r1/_R1.fastq/_R2.fastq}
                r2=\${r2/_1.fastq/_2.fastq}
                if [[ "\$r1" == *".gz" ]]; then r2=\${r2}.gz; fi

                bowtie2 -x ./\${name}_db -1 \$r1 -2 \$r2 -S tmp.sam -p ${task.cpus}
                samtools view -b -F 12 tmp.sam | samtools fastq -1 tmp_1.fq -2 tmp_2.fq -0 /dev/null -s /dev/null -
                cat tmp_1.fq >> \${name}_mapped_1.fq
                cat tmp_2.fq >> \${name}_mapped_2.fq
                rm -f tmp.sam tmp_1.fq tmp_2.fq
            elif [[ "\$fq" != *_R2.fastq* && "\$fq" != *_2.fastq* ]]; then
                r1=\$fq
                bowtie2 -x ./\${name}_db -U \$r1 -S tmp.sam -p ${task.cpus}
                samtools view -b -F 12 tmp.sam | samtools fastq -0 /dev/null -s /dev/null - > tmp_single.fq
                cat tmp_single.fq >> \${name}_mapped_single.fq
                rm -f tmp.sam tmp_single.fq
            fi
        done
        if [ -s "\${name}_mapped_1.fq" ]; then
            echo "Running Paired-End assembly for \$name"
            spades.py --sc --careful -1 ./\${name}_mapped_1.fq -2 ./\${name}_mapped_2.fq -o reassembled_bins/\${name}_out --threads ${task.cpus}
        elif [ -s "\${name}_mapped_single.fq" ]; then
            echo "Running Single-End assembly for \$name"
            spades.py --sc --careful -s ./\${name}_mapped_single.fq -o reassembled_bins/\${name}_out --threads ${task.cpus}
        else
            echo "Warning: No reads extracted for target \$name, skipping assembly."
        fi
        rm -rf ./target.fa ./\${name}_db ./\${name}_mapped_1.fq ./\${name}_mapped_2.fq ./\${name}_mapped_single.fq
    done
    """
}
