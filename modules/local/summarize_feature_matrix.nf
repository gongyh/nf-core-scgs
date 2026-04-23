process SUMMARIZE_FEATURE_MATRIX {
    tag "summarize_all"
    label 'process_low'
    conda "bioconda::hamronization=1.1.9"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hamronization:1.1.9--pyhdfd78af_0'
        : 'biocontainers/hamronization:1.1.9--pyhdfd78af_0'}"


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
    python3 - <<'EOF'
import csv
data_depth = {}
header_depth = []
with open('${depth}', 'r') as f:
    reader = csv.reader(f, delimiter='\\t')
    header_depth = next(reader)
    for row in reader:
        if row:
            data_depth[row[0]] = row[1:]
data_kmer = {}
header_kmer = []
with open('${k4_csv}', 'r') as f:
    reader = csv.reader(f, delimiter=',')
    header_kmer = next(reader)
    for row in reader:
        if row:
            data_kmer[row[0]] = row[1:]

common_ids = sorted(set(data_depth.keys()) & set(data_kmer.keys()))

with open('final_feature_matrix.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([header_depth[0]] + header_depth[1:] + header_kmer[1:])
    for cid in common_ids:
        writer.writerow([cid] + data_depth[cid] + data_kmer[cid])
EOF

cat <<EOF > versions.yml
"${task.process}":
    python: \$(python3 --version | sed 's/Python //')
EOF
    """
}
