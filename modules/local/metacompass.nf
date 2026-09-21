nextflow.enable.types = true

process METACOMPASS {
    tag "${meta.id}"
    label 'process_high'

    conda "scgs::metacompass=1.12=1.12--h9948957_6 bioconda::seqkit=2.10.0 bioconda::seqtk=1.4"
    container "scgs/mulled-v2-0e7fe6bd3265990ffcdf96496fe08dc5aa55fd24:62c3137bd1d05677122f8069cb3981ac4e60651e-6"

    input:
    tuple(meta: Map, reads: List<Path>, refs_fna: List<Path>)

    output:
    record(meta: meta, contig: file("${prefix}_*.metacompass.ctg.fa"), assembly: file("${prefix}_*.metacompass_out", type: 'dir'))
    topic:
    file('versions.yml') >> 'versions'

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
