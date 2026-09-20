nextflow.enable.types = true

process DOUBLETD {
    label 'process_medium'

    conda "pypy=7.3.3 python=3.6.12 numpy=1.19.5 scipy=1.5.3 pandas=1.1.5 pyvcf=0.6.8"
    container "community.wave.seqera.io/library/pypy_python_numpy_scipy_pruned:a4047355723aa61e"

    input:
    vcf: Path

    output:
    record(prediction: file('prediction.tsv'))
    topic:
    file('versions.yml') >> 'versions'

    script:
    """
    vcf_csv.py
    doubletd_haploid.py --inputAlternate AD.csv --inputTotal DP.csv --delta 0.2 --beta 0.05 -o prediction.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        doubletD: 0.1.0
    END_VERSIONS
    """
}
