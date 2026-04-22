process SAMTOOLS_COVERAGE_COMBINED {
    tag "all_samples"
    label 'process_medium'
    conda "bioconda::samtools=1.17 conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00c71ee_0' :
        'docker.io/biocontainers/samtools:1.17--h00c71ee_0' }"
    input:
    path bams
    path fasta
    path fai
    output:
    path "abundance_matrix.tsv", emit: matrix
    path "versions.yml"        , emit: versions

    script:
    """
    for bam_file in ${bams}; do
        sample_name=\$(basename \${bam_file} .bam)
        samtools coverage --reference ${fasta} -o \${sample_name}.cov \${bam_file}
    done
        python3 -c "
import csv
import glob
cov_files = glob.glob('*.cov')
samples = []
data = {}  # contig_id -> {sample: meandepth}
for f in cov_files:
    sample = f.replace('.cov', '')
    samples.append(sample)
    with open(f, 'r') as inf:
        reader = csv.reader(inf, delimiter='\\t')
        header = next(reader)
        for row in reader:
            contig = row[0]
            meandepth = float(row[6])
            if contig not in data:
                data[contig] = {}
            data[contig][sample] = meandepth
with open('abundance_matrix.tsv', 'w', newline='') as outf:
    writer = csv.writer(outf, delimiter='\\t')
    writer.writerow(['contig_id'] + samples)
    for contig in sorted(data.keys()):
        row = [contig] + [str(data[contig].get(s, '0')) for s in samples]
        writer.writerow(row)
"

    cat <<EOF > versions.yml
    "${task.process}":
        samtools: \$(samtools version | sed '1!d;s/.* //')
        python: \$(python3 --version | sed 's/Python //')
    EOF
    """
}
