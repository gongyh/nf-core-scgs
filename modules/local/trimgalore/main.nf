process TRIMGALORE {
    tag "$name"
    label 'process_high'

    publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

    conda "bioconda::trim-galore=0.6.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.7--hdfd78af_0' :
        'biocontainers/trim-galore:0.6.7--hdfd78af_0' }"

    input:
    tuple val(name), file(reads)
    val(single_end)

    output:
    path("*.fq.gz"),                      emit: trimmed_reads
    path("*trimming_report.txt"),         emit: results
    path("*_fastqc.{zip,html}"),          emit: reports

    script:
    def c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    def c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    def tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    def tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    if (single_end) {
        """
        trim_galore --trim-n --max_n 0 --fastqc --gzip --fastqc_args \"--threads ${task.cpus}\" --cores 4 $c_r1 $tpc_r1 $reads
        """
    } else {
        """
        trim_galore --paired --trim-n --max_n 0 --fastqc --gzip --fastqc_args \"--threads ${task.cpus}\" --cores 4 $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    }
}
