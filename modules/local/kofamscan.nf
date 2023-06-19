process KOFAMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kofamscan=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kofamscan:1.3.0--hdfd78af_2':
        'biocontainers/kofamscan:1.3.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(faa)
    path profile
    path ko_list

    output:
    tuple val(meta), path("${prefix}_KOs_*.txt"), emit: txt

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -o ${prefix}_KOs_detail.txt ${faa}
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -r -f mapper -o ${prefix}_KOs_mapper.txt ${faa}
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -r -f mapper-one-line -o ${prefix}_KOs_mapper2.txt ${faa}
    kofam_postprocess.py /opt/nf-core-scgs/assets/ko_KO.txt ${prefix}_KOs_mapper.txt > ${prefix}_KOs_ko.txt
    """
}
