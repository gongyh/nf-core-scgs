process CANU {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::canu=2.2 bioconda::flye=2.9.2 bioconda::perl-bioperl=1.7.8 bioconda::perl-io-all=0.87"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ef18b85941472064953d94d6112dc05e36194472:90a061fb09200847c4637e3480aa75a467fb2786-0':
        'scgs/mulled-v2-ef18b85941472064953d94d6112dc05e36194472:90a061fb09200847c4637e3480aa75a467fb2786-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.ctg200.fasta"), emit: ctg200
    tuple val(meta), path("${prefix}.ctgs.fasta")  , emit: ctg
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def mode = params.bulk ? "bulk" : "mda"
    """
    if [ \"${mode}\" == \"bulk\" ]; then
        #canu -d ${prefix}.spades_out -p ${prefix} genomeSize=4m useGrid=false maxThreads=${task.cpus} maxMemory=${task.memory.toGiga()}g -nanopore ${reads[0]}
        flye --nano-raw ${reads[0]} --out-dir ${prefix}.spades_out --threads ${task.cpus} --scaffold
    else
        #canu -d ${prefix}.spades_out -p ${prefix} genomeSize=4m corOutCoverage=999 corMhapSensitivity=high useGrid=false maxThreads=${task.cpus} maxMemory=${task.memory.toGiga()}g -nanopore ${reads[0]}
        flye --nano-raw ${reads[0]} --out-dir ${prefix}.spades_out --threads ${task.cpus} --scaffold --meta
    fi
    #ln -s ${prefix}.spades_out/${prefix}.contigs.fasta ${prefix}.contigs.fasta # for canu
    cut -f1,2,3 ${prefix}.spades_out/assembly_info.txt | awk -F'\t' 'NR>1{print \$1"\t"\$1"_length_"\$2"_cov_"\$3}' > flyeID_spadesID.txt
    fasta_tool --swap_ids flyeID_spadesID.txt ${prefix}.spades_out/assembly.fasta > ${prefix}.contigs.fasta
    ##ln -s ${prefix}.spades_out/assembly.fasta ${prefix}.contigs.fasta # for flye
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/ len=.*\$//g' > ${prefix}.ctgs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        canu: \$(echo \$(canu -version 2>&1) | sed 's/^.*canu //; s/Using.*\$//')
    END_VERSIONS
    """
}
