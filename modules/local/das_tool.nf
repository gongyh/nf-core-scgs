process DAS_TOOL {
    tag "das_tool"
    label 'process_medium'
    conda "bioconda::das_tool=1.1.2"
    container 'community.wave.seqera.io/library/das_tool:1.1.2--0fc15370c91e86b2'

    input:
    path assembly
    val  raw_info


    output:
    path "das_tool_bins", emit: bins
    path "versions.yml", emit: versions

    script:
    def my_labels = raw_info.findAll { it instanceof String }.join(',')
    def my_files  = raw_info.findAll { it instanceof nextflow.processor.TaskPath || it instanceof java.nio.file.Path }.join(',')
    def args = task.ext.args ?: ''
    """
    DAS_Tool -i ${my_files} \\
        -c ${assembly} \\
        -l ${my_labels} \\
        -o das_tool_result \\
        --search_engine diamond \\
        --threads ${task.cpus} \\
        ${args}

    mkdir -p das_tool_bins
    if [ -d das_tool_result_DASTool_bins ]; then
        cp das_tool_result_DASTool_bins/*.fa das_tool_bins/ 2>/dev/null || true
    elif [ -d das_tool_result_bins ]; then
        cp das_tool_result_bins/*.fa das_tool_bins/ 2>/dev/null || true
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        das_tool: \$(DAS_Tool --version 2>&1 | head -1)
    END_VERSIONS
    """
}

