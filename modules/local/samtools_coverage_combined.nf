process SAMTOOLS_COVERAGE_COMBINED {
    tag "all_samples"
    label 'process_medium'
    publishDir "${params.outdir}/coverage_matrix", mode: 'copy'
    conda "bioconda::samtools=1.17 conda-forge::pandas=1.5.3 conda-forge::python=3.11"

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

    python3 -c "
    import pandas as pd
    import glob
    import os

    all_files = glob.glob('*.cov')
    combined_data = {}

    for f in all_files:
        sample = f.replace('.cov', '')
        df = pd.read_csv(f, sep='\\t')
        combined_data[sample] = df.set_index('#rname')['meandepth']

    matrix = pd.DataFrame(combined_data)
    matrix.index.name = 'contig_id'
    matrix.to_csv('abundance_matrix.tsv', sep='\\t')
    "

    cat <<EOF > versions.yml
    "${task.process}":
        samtools: \$(samtools version | sed '1!d;s/.* //')
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)')
    EOF
    """
}
