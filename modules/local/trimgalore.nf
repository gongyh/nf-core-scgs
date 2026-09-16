nextflow.enable.types = true

process TRIMGALORE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::trim-galore=0.6.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.7--hdfd78af_0' :
        'biocontainers/trim-galore:0.6.7--hdfd78af_0' }"

    input:
    tuple(meta: Map, reads: List<Path>)

    output:
    record(meta: meta, reads: file('*{3prime,5prime,trimmed,val}{,_1,_2}.fq.gz'), log: file('*trimming_report.txt', optional: true), zip: file('*.zip', optional: true), versions: file('versions.yml'))

    script:
    def c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    def c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    def tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    def tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    if (meta.single_end) {
        """
        trim_galore --trim-n --max_n 0 --fastqc --gzip --fastqc_args \"--threads ${task.cpus}\" --cores 4 ${c_r1} ${tpc_r1} ${reads}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    } else {
        """
        trim_galore --paired --trim-n --max_n 0 --fastqc --gzip --fastqc_args \"--threads ${task.cpus}\" --cores 4 ${c_r1} ${c_r2} ${tpc_r1} ${tpc_r2} ${reads}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
}
