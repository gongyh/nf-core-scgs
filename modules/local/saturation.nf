process SATURATION {
    tag "${meta.id}"

    conda "fastp=0.20.1,mccortex=1.0,r-base=4.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-be9f5ba45af47384627efd3880c6586bd81bf92a:2d818d847653deccb4eefa54656b340b6eccd5b7-0' :
        'scgs/mulled-v2-be9f5ba45af47384627efd3880c6586bd81bf92a:2d818d847653deccb4eefa54656b340b6eccd5b7-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    path("${prefix}_kmer.pdf")
    path("${prefix}_cov31_*.csv")

    when:
    params.saturation

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
    """
    fastp -i ${reads[0]} -A -G -Q -L -s 10 -d 0 -o ${prefix}_split.fq.gz
    for i in {1..10}; do
    mccortex31 build --kmer 31 --sample \$i -t ${task.cpus} -Q 20 -m 8G \
        --seq \${i}.${prefix}_split.fq.gz \${i}.k31.ctx
    if [ \$i == 1 ]; then
        mccortex31 clean -t ${task.cpus} -m 8G -U10 -T16 -f -o null -C ${prefix}_cov31_p\${i}.csv 0:\${i}.k31.ctx
        cp -f 1.k31.ctx tmp_clean31.ctx
    else
        mccortex31 join -m 8G --out merged_clean31.ctx 0:\${i}.k31.ctx 0:tmp_clean31.ctx
        mccortex31 clean -t ${task.cpus} -m 8G -U10 -T16 -f -o null -C ${prefix}_cov31_p\${i}.csv 0:merged_clean31.ctx
        mv -f merged_clean31.ctx tmp_clean31.ctx
    fi
    done
    KmerDensity.R \$PWD ${prefix}
    """
    } else {
    """
    fastp -i ${reads[0]} -I ${reads[1]} -A -G -Q -L -s 10 -d 0 -o ${prefix}_split_R1.fq.gz -O ${prefix}_split_R2.fq.gz
    for i in {1..10}; do
    mccortex31 build --kmer 31 --sample \$i -t ${task.cpus} -Q 20 -m 8G \
        --seq2 \${i}.${prefix}_split_R1.fq.gz:\${i}.${prefix}_split_R2.fq.gz \${i}.k31.ctx
    if [ \$i == 1 ]; then
        mccortex31 clean -t ${task.cpus} -m 8G -U10 -T16 -f -o null -C ${prefix}_cov31_p\${i}.csv 0:\${i}.k31.ctx
        cp -f 1.k31.ctx tmp_clean31.ctx
    else
        mccortex31 join -m 8G --out merged_clean31.ctx 0:\${i}.k31.ctx 0:tmp_clean31.ctx
        mccortex31 clean -t ${task.cpus} -m 8G -U10 -T16 -f -o null -C ${prefix}_cov31_p\${i}.csv 0:merged_clean31.ctx
        mv -f merged_clean31.ctx tmp_clean31.ctx
    fi
    done
    KmerDensity.R \$PWD ${prefix}
    """
    }

}
