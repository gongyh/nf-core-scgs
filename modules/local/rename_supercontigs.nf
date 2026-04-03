process RENAME_SUPERCONTIGS {
    tag "rename"
    label 'process_low'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path contigs

    output:
    path "super_contigs.fasta", emit: super_contigs
    path "super_contigs_name.txt", emit: name_map
    path "versions.yml", emit: versions

    script:
    """
    python3 <<EOF
import sys
with open("${contigs}", "r") as fin, open("super_contigs.fasta", "w") as fout, open("super_contigs_name.txt", "w") as fmap:
    counter = 1
    for line in fin:
        if line.startswith(">"):
            old_name = line[1:].strip()
            new_name = f"SuperContig_{counter}"
            fout.write(f">{new_name}\\n")
            fmap.write(f"{new_name}\\t{old_name}\\n")
            counter += 1
        else:
            fout.write(line)
EOF
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        script: "rename_supercontigs.nf"
    END_VERSIONS
    """
}
