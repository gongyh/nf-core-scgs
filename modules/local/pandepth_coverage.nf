process PANDEPTH_COVERAGE {
    tag "${bam.getBaseName()}"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${bam.getBaseName()}.depth"), emit: depth
    path "versions.yml"                                 , emit: versions

    script:
    def pandepth_bin = "${projectDir}/bin/pandepth"
    def threads = task.cpus ?: 1
    def sample_name = bam.getBaseName()
    """
    if [ ! -f "${bam}.bai" ]; then
        samtools index "${bam}"
    fi
    ${pandepth_bin} -i "${bam}" -t ${threads} -r "${fasta}" -o "${sample_name}"
    zcat "${sample_name}.chr.stat.gz" 2>/dev/null | awk 'NR>1 {print \$1"\\t"\$5}' | sort -k1,1 > "${sample_name}.depth"
    rm -f "${sample_name}.chr.stat.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandepth: \$(${pandepth_bin} -h 2>&1 | head -1)
    END_VERSIONS
    """
}
