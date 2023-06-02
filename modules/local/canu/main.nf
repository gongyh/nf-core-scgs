process CANU {
    tag "${prefix}"
    label 'process_high'
    publishDir "${params.outdir}/spades", mode: 'copy'

    conda "bioconda::canu=2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/canu:2.2--ha47f30e_0':
        'biocontainers/canu:2.2--ha47f30e_0' }"

    input:
    path clean_reads

    output:
    path("${prefix}.ctg200.fasta"),                     emit: ctg200
    path("${prefix}.ctgs.fasta"),                       emit: ctg

    when:
    params.ass

    script:
    prefix = clean_reads[0].toString() - ~/(_trimmed)?(_norm)?(_combined)?(\.R1)?(_1)?(_R1)?(\.1_val_1)?(_1_val_1)?(_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?(\.bz2)?$/
    def R1 = clean_reads[0].toString()
    def mode = params.bulk ? "bulk" : "mda"
    """
    if [ \"${mode}\" == \"bulk\" ]; then
    #canu -d ${prefix}.spades_out -p ${prefix} genomeSize=4m useGrid=false maxThreads=${task.cpus} maxMemory=${task.memory.toGiga()}g -nanopore $R1
    flye --nano-raw $R1 --out-dir ${prefix}.spades_out --threads ${task.cpus} --scaffold
    else
    #canu -d ${prefix}.spades_out -p ${prefix} genomeSize=4m corOutCoverage=999 corMhapSensitivity=high useGrid=false maxThreads=${task.cpus} maxMemory=${task.memory.toGiga()}g -nanopore $R1
    flye --nano-raw $R1 --out-dir ${prefix}.spades_out --threads ${task.cpus} --scaffold --meta
    fi
    #ln -s ${prefix}.spades_out/${prefix}.contigs.fasta ${prefix}.contigs.fasta # for canu
    cut -f1,2,3 ${prefix}.spades_out/assembly_info.txt | awk -F'\t' 'NR>1{print \$1"\t"\$1"_length_"\$2"_cov_"\$3}' > flyeID_spadesID.txt
    fasta_tool --swap_ids flyeID_spadesID.txt ${prefix}.spades_out/assembly.fasta > ${prefix}.contigs.fasta
    ##ln -s ${prefix}.spades_out/assembly.fasta ${prefix}.contigs.fasta # for flye
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/ len=.*\$//g' > ${prefix}.ctgs.fasta
    """
}
