process UMAP {
    label 'process_medium'

    conda "conda-forge::scanpy=1.11.2 conda-forge::pyyaml=6.0.2"
    container "scgs/mulled-v2-9109a57576476a9373a70c6f48e5d8d64c8d6c77:1da6c154f4395390d27ed0224ea5055605acd644-0"

    input:
    path("tda/*")

    output:
    path "umap.h5ad"
    path "umap.pkl"
    path "umap.pdf"
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('umap.py')

    stub:
    """
    touch "umap.h5ad"
    touch "umap.pkl"
    touch "umap.pdf"
    touch "versions.yml"
    """
}
