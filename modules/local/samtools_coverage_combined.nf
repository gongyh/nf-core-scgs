process SAMTOOLS_COVERAGE_COMBINED {
    tag "all_samples"
    label 'process_medium'
    conda "bioconda::samtools=1.17 conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00c71ee_0' :
        'docker.io/biocontainers/samtools:1.17--h00c71ee_0' }"
    input:
    path bams
    path fasta
    path fai
    output:
    path "abundance_matrix.tsv", emit: matrix
    path "versions.yml"        , emit: versions

    script:
    """
    for bam_file in ${bams}; do
        sample_name=\$(basename \${bam_file} .bam)
        samtools coverage --reference ${fasta} -o \${sample_name}.cov \${bam_file}
    done
    for cov in *.cov; do
        sample=\${cov%.cov}
        awk '!/^#/ {print \$1"\t"\$7}' \$cov | sort -k1,1 > \${sample}.depth
    done
    samples=(\$(ls *.depth | sed 's/.depth//'))
    cut -f1 *.depth | sort -u > all_contigs.tmp
    for sample in \${samples[*]}; do
        join -a1 -e0 -o '2.2' -t \$'\t' all_contigs.tmp \${sample}.depth > \${sample}.depth_col
    done
    paste all_contigs.tmp \$(for s in \${samples[*]}; do echo \${s}.depth_col; done) > abundance_matrix.tsv
    header="contig_id"
    for sample in \${samples[*]}; do
        header="\${header}\t\${sample}"
    done
    (echo -e "\${header}" && cat abundance_matrix.tsv) > abundance_matrix.tsv.tmp && mv abundance_matrix.tsv.tmp abundance_matrix.tsv
    rm -f *.depth *.depth_col all_contigs.tmp
    cat <<EOF > versions.yml
    "${task.process}":
        samtools: \$(samtools version | sed '1!d;s/.* //')
    EOF
    """
}
