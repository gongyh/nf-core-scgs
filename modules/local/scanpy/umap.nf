process UMAP {
    label 'process_medium'

    conda "conda-forge::scanpy=1.11.2 conda-forge::pyyaml=6.0.2 conda-forge::python-igraph=0.11.9 conda-forge::leidenalg=0.10.2"
    container "scgs/mulled-v2-b72f30682a2b6401a020a481725ded0634ad5f6c:2ebc80117a8ef44f0d77a014d339c39059767495-0"

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
