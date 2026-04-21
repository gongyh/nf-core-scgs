process KMER_COUNT {
    tag "$meta.id - k$kmer"
    label 'process_low'
    publishDir "${params.outdir}/kmer", mode: 'copy'
    conda "conda-forge::pandas=1.5.3 conda-forge::biopython=1.81 conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python-pandas-biopython:1.0' :
        'evolbioinfo/python-pandas-biopython:1.0' }"

    input:
    tuple val(meta), path(fasta)
    val kmer

    output:
    tuple val(meta), val(kmer), path("${meta.id}_k${kmer}.csv"), emit: csv
    path "versions.yml"                                         , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python - <<EOF
    import pandas as pd
    from Bio import SeqIO
    from itertools import product
    import collections
    import sys

    def get_kmers(k):
        return [''.join(p) for p in product('ACGT', repeat=k)]

    def count_kmers(fasta_file, k):
        kmers = get_kmers(k)
        results = []

        for record in SeqIO.parse(fasta_file, "fasta"):
            name = record.id
            sequence = str(record.seq).upper()

            counts = collections.Counter()
            for i in range(len(sequence) - k + 1):
                kmer_seq = sequence[i:i+k]
                if 'N' not in kmer_seq and len(kmer_seq) == k:
                    counts[kmer_seq] += 1

            total = sum(counts.values()) if sum(counts.values()) > 0 else 1
            row = {'contig_id': name}
            for k_str in kmers:
                row[k_str] = counts[k_str] / total
            results.append(row)

        if not results:
            df = pd.DataFrame(columns=['contig_id'] + kmers)
        else:
            df = pd.DataFrame(results)

        df.to_csv('${prefix}_k${kmer}.csv', index=False)

    try:
        count_kmers('$fasta', int('$kmer'))
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
