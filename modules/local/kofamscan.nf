process KOFAMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kofamscan=1.3.0 conda-forge::python=3.6.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ef3cc10895f39bdde312c5e796de361bc231bb29:f6fe8bf9968d952a4cb8cdb90f165e49c40688b8-0':
        'scgs/mulled-v2-ef3cc10895f39bdde312c5e796de361bc231bb29:f6fe8bf9968d952a4cb8cdb90f165e49c40688b8-0' }"

    input:
    tuple val(meta), path(faa)
    path profile
    path ko_list

    output:
    tuple val(meta), path("${prefix}_KOs_*.txt"), emit: txt
    path "versions.yml"                         , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -o ${prefix}_KOs_detail.txt ${faa}
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -r -f mapper -o ${prefix}_KOs_mapper.txt ${faa}
    exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -r -f mapper-one-line -o ${prefix}_KOs_mapper2.txt ${faa}
    kofam_postprocess.py \$(echo \$(which ko_KO.txt)) ${prefix}_KOs_mapper.txt > ${prefix}_KOs_ko.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        koramscan: \$(echo \$(exec_annotation -v 2>&1) | sed 's/^.*exec_annotation //; s/Using.*\$//')
    END_VERSIONS
    """
}
