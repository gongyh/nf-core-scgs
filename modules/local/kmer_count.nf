process KMER_COUNT {
    tag "$meta.id - k$kmer"
    label 'process_low'
    conda "conda-forge::opentsne=1.0.0 conda-forge::h5py=3.9.0 conda-forge::numpy=1.25.0 conda-forge::pandas=2.0.2 bioconda::kpal=2.1.1 bioconda::perl-bioperl=1.7.8"
    container "scgs/mulled-v2-8905087433117c98a93e379c07447431e85bdd71:5402918794aa21f8f7e4b46973655d86142c9ffb-0"

    input:
    tuple val(meta), path(fasta)
    val kmer

    output:
    tuple val(meta), val(kmer), path("${meta.id}_k${kmer}.csv"), emit: csv
    path "versions.yml"                                         , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def k = kmer
    """
    python3 - <<EOF
    import pandas as pd
    from Bio import SeqIO
    import kpal.klib as klib
    from itertools import product
    import numpy as np
    import sys

    def run_kpal_matrix(infile, k, outfile):
        kmers = [''.join(p) for p in product('ACGT', repeat=k)]
        results = []
        contig_ids = []
        for record in SeqIO.parse(infile, "fasta"):
            seq = str(record.seq).upper()
            if len(seq) < k:
                continue
            try:
                counts = klib.count(seq, k)
                total = counts.sum()
                freqs = counts / total if total > 0 else counts
                results.append(freqs)
                contig_ids.append(record.id)
            except Exception as e:
                print(f"Warning: Could not process {record.id}: {e}", file=sys.stderr)
        if results:
            df = pd.DataFrame(results, columns=kmers)
            df.insert(0, 'contig_id', contig_ids)
            df.to_csv(outfile, index=False)
        else:
            pd.DataFrame(columns=['contig_id'] + kmers).to_csv(outfile, index=False)

    run_kpal_matrix("${fasta}", int("${k}"), "${prefix}_k${k}.csv")
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        kpal: \$(kpal --version 2>&1 | sed 's/kpal //')
    END_VERSIONS
    """
}
