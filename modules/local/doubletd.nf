process DOUBLETD {
    label 'process_medium'

    conda "conda-forge::numpy=1.25.0 bioconda::scipy=1.10.1 bioconda::pysam=0.21.0 bioconda::pandas=2.0.2 bioconda::pyvcf=0.6.8 conda-forge::python=3.9.16"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-cabd179211b935b76ccbd0bcec69e95eab67a374:aff7fa7836f60ba389a388528d6cfaf626926614-0' :
        'scgs/mulled-v2-cabd179211b935b76ccbd0bcec69e95eab67a374:aff7fa7836f60ba389a388528d6cfaf626926614-0' }"

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
