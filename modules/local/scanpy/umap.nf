nextflow.enable.types = true

process UMAP {
    label 'process_medium'

    conda "conda-forge::scanpy=1.11.2 conda-forge::pyyaml=6.0.2 conda-forge::python-igraph=0.11.9 conda-forge::leidenalg=0.10.2 conda-forge::plotly=6.2.0"
    container "scgs/mulled-v2-bebad6fb9c0a642cb203291e2b9969552cec05d6:955d7191e655a067018f09dfc80d57c23afb23c9-0"

    input:
    tda: Bag<Path>

    output:
    record(h5ad: file('umap.h5ad'), pkl: file('umap.pkl'), pdf: file('umap.pdf'), html: file('umap.html'), versions: file('versions.yml'))
    topic:
    file('versions.yml') >> 'local_versions'

    script:
    template('umap.py')

    stub:
    """
    touch "umap.h5ad"
    touch "umap.pkl"
    touch "umap.pdf"
    touch "umap.html"
    touch "versions.yml"
    """
}
