process KOFAMSCAN {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(faa)
    path profile
    path ko_list

    output:
    path("${prefix}_KOs_*.txt"), emit: txt

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -o ${prefix}_KOs_detail.txt ${faa}
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -r -f mapper -o ${prefix}_KOs_mapper.txt ${faa}
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -r -f mapper-one-line -o ${prefix}_KOs_mapper2.txt ${faa}
    kofam_postprocess.py /opt/nf-core-scgs/assets/ko_KO.txt ${prefix}_KOs_mapper.txt > ${prefix}_KOs_ko.txt
    """
}
