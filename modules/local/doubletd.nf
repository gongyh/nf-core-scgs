process DOUBLETD {
    label 'process_medium'

    conda "pypy=7.3.3 python=3.6.12 numpy=1.19.5 scipy=1.5.3 pysam=0.15.4 pandas=1.1.5 pyvcf=0.6.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b7b08ecbb8cbfa2e328c1c68384966d9ba24b2c8:35b29e00884f8c94506f5174a2a2eb3ea2f07a20-0' :
        'scgs/mulled-v2-b7b08ecbb8cbfa2e328c1c68384966d9ba24b2c8:35b29e00884f8c94506f5174a2a2eb3ea2f07a20-0' }"

    input:
    path("monovar.vcf")

    output:
    path("prediction.tsv"), emit: prediction
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

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
