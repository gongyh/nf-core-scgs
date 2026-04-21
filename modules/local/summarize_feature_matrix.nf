process SUMMARIZE_FEATURE_MATRIX {
    tag "summarize_all"
    label 'process_low'
    publishDir "${params.outdir}/feature_matrix", mode: 'copy'
    conda "bioconda::samtools=1.17 conda-forge::pandas=1.5.3 conda-forge::python=3.11"
    container "community.wave.seqera.io/library/samtools_pandas:bc6974910398686e"
    input:
    tuple val(meta), path(fasta)
    path depth
    tuple val(meta2), val(kmer_size), path(k4_csv)
    tuple val(meta3), path(gff)

    output:
    path "final_feature_matrix.csv", emit: matrix
    path "versions.yml"            , emit: versions

    script:
    """
    python -c "
    import pandas as pd

    df_depth = pd.read_csv('${depth}', sep='\\t', index_col=0)

    df_kmer = pd.read_csv('${k4_csv}', index_col=0)

    result = pd.concat([df_depth, df_kmer], axis=1)


    result.to_csv('final_feature_matrix.csv')
    "

    cat <<EOF > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)')
    EOF
    """
}
