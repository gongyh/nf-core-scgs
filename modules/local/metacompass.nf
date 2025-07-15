process METACOMPASS {
    tag "${meta.id}"
    label 'process_high'

    conda "scgs::metacompass=1.12=1.12--h9948957_5 bioconda::seqkit=2.10.0 bioconda::seqtk=1.4"
    container "scgs/mulled-v2-0e7fe6bd3265990ffcdf96496fe08dc5aa55fd24:62c3137bd1d05677122f8069cb3981ac4e60651e-1"

    input:
    tuple val(meta), path(reads)
    path(refs_fna)

    output:
    tuple val(meta), path("${prefix}_*.metacompass.ctg.fa")                , emit: contig
    tuple val(meta), path("${prefix}_*.metacompass_out")                   , emit: assembly
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def rcl = meta.single_end ? "-U ${reads[0]}" : "-P ${reads[0]},${reads[1]}"
    """
    refs_fna=(${refs_fna})
    for ref_fna in \${refs_fna[*]}; do
        go_metacompass.py -r \${ref_fna} ${rcl} -m 1 -g 100 -e ${prefix} -t ${task.cpus} -o ${prefix}_\${ref_fna}.metacompass_out
        cp ${prefix}_\${ref_fna}.metacompass_out/metacompass.final.ctg.fa ${prefix}_\${ref_fna}.metacompass.ctg.fa
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaCompass: 1.12
    END_VERSIONS
    """
}
