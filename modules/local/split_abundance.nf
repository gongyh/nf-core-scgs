process SPLIT_ABUNDANCE {
    tag "split"
    label 'process_low'

    conda "bioconda::maxbin2=2.2.7"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/maxbin2:2.2.7--he1b5a44_2' :
        'quay.io/biocontainers/maxbin2:2.2.7--he1b5a44_2' }"
    input:
    path coverage_tsv

    output:
    path "*.abund", emit: abund_files

    script:
    """
    awk 'NR==1 {
        for(i=2;i<=NF;i++) {
            col[i]=\$i;
            out=col[i]".abund";
            print "contig_id\\t"col[i] > out
        }
    }
    NR>1 {
        for(i=2;i<=NF;i++) {
            out=col[i]".abund";
            print \$1"\\t"\$i >> out
        }
    }' ${coverage_tsv}
    """
}
