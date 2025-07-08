process UMAP {
    label 'process_medium'

    conda "bioconda::scanpy-scripts=1.9.301 conda-forge::python=3.9.23 conda-forge::loompy=3.0.6 conda-forge::pyyaml=6.0.2 conda-forge::r-ggpubr=0.6.1 bioconda::r-sceasy=0.0.7"
    container "scgs/mulled-v2-fe9371f6be95d197dedf7c0a65e9e322526829e0:8d6c60157f466e2d016b24eb850b151a57f390ab-0"

    input:
    path("tda/*", arity: '3..*')

    output:
    path "umap.h5ad"
    path "umap.pkl"
    path "umap.pdf"
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('umap.py')
    """
    scanpy-cli plot embed --projection 2d --color sample_genus --title UMAP umap.h5ad umap.pdf
    """

    stub:
    """
    touch "umap.h5ad"
    touch "umap.pkl"
    touch "versions.yml"
    """
}
