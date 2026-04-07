// 完整流程：输入 → 合并 → 质量修剪 → 子样本组装 → 联合组装 → 比对 → 覆盖度 → 共现分箱
// 最终输出：super_contigs.fasta，共现分箱结果（clusters.tsv, bins/）

// 显示帮助信息（当用户输入 --help 时）
def helpMessage() {
    // 打印 nf-core 风格的头部信息
    log.info nfcoreHeader()
    // 打印使用说明
    log.info"""

    Usage:

    The typical command for running the minimeta pipeline is as follows:

    nextflow run gongyh/nf-core-scgs --reads '*_R{1,2}.fastq.gz' --minimeta -profile docker

    Mandatory arguments:
    --reads                       Path to input data (must be surrounded with quotes)
    -profile                      Configuration profile to use. Can use multiple (comma separated). Available: conda, docker, singularity, awsbatch, test and more.

    Options:
    --single_end                  Specifies that the input is single end reads
    --notrim                      Specifying --notrim will skip the adapter trimming step.
    --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.

    Trimming options:
    --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1
    --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2
    --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1
    --three_prime_clip_r2 [int]   Instructs Trim Galore to remove bp from the 3' end of read 2

    Output options:
    --outdir                      The output directory where the results will be saved
    --email                       Set this parameter to your e-mail address to get a summary e-mail
    --maxMultiqcEmailFileSize     Threshold size for MultiQC report to be attached in notification email (Default: 25MB)

    AWSBatch options:
    --awsqueue                    The AWSBatch JobQueue
    --awsregion                   The AWS Region
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 * 初始化配置变量
 */

// 如果用户指定了 --help，则显示帮助信息并退出
if (params.help){
    helpMessage()
    exit 0
}

// 设置默认参数值
params.single_end = false          // 默认双端测序
params.notrim = false              // 默认进行质量修剪
params.saveTrimmed = false         // 默认不单独保存修剪后的文件

// 当前运行的名称（用于日志）
custom_runName = workflow.runName
// 单端标志（局部变量）
single_end = params.single_end

// AWSBatch 相关配置检查（如果使用 awsbatch 配置文件）
if(workflow.profile == 'awsbatch') {
    // 检查必须的 aws 参数
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // 检查工作目录和输出目录是否在 S3 上
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// 创建 MultiQC 配置文件的通道（如果用户提供）
ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)
// 空的通道，用于扩展
ch_multiqc_custom_config = Channel.empty()
ch_multiqc_logo = Channel.empty()
// 输出文档的通道
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

// 自定义修剪选项（默认均为 0，表示不额外裁剪）
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

/*
 * Create a channel for input read files
 * 创建输入读取文件的通道
 */
if(params.readPaths){
    // 如果用户提供了 readPaths（直接指定文件路径列表）
    if(single_end){
        // 单端模式
        read_files_fastqc = read_files_trimming =
        Channel.from(params.readPaths, checkIfExists: true)
            .map { row -> def meta=[:];
                    meta.id = row[0];                       // 样本ID
                    meta.single_end = single_end;           // 单端标志
                    [meta, [file(row[1][0]), file(row[1][1])]]}  // 输出 [meta, [R1, R2]]
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
    } else {
        // 双端模式
        read_files_fastqc = read_files_trimming =
        Channel.from(params.readPaths)
            .map { row -> def meta=[:];
                    meta.id = row[0];
                    meta.single_end = single_end;
                    [meta, [file(row[1][0]), file(row[1][1])]]}
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
    }
} else {
    // 使用 --reads 参数自动匹配文件对
    if (single_end) {
        // 单端模式：每个文件单独处理
        read_files_fastqc = read_files_trimming =
        Channel.fromFilePairs(params.reads, size:1, checkIfExists: true)
            .map { it ->
                def meta = [:];
                meta.id = it[0].replaceFirst(~/\.[^\.]+$/, '');   // 去除扩展名作为ID
                meta.single_end = single_end;
                [meta, [file(it[1][0])]]}
    } else {
        // 双端模式：自动匹配 R1/R2 对（默认模式）
        read_files_fastqc = read_files_trimming =
        Channel.fromFilePairs(params.reads, size:2, checkIfExists: true)
            .map { it ->
                def meta = [:];
                meta.id = it[0].replaceFirst(~/\.[^\.]+$/, '');
                meta.single_end = single_end;
                [meta, [file(it[1][0]), file(it[1][1])]]}
    }
}

// 打印运行摘要信息
log.info nfcoreHeader()
def summary = [:]                          // 创建摘要映射
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']            = params.reads
summary['Data Type']        = single_end ? 'Single-End' : 'Paired-End'
summary['Workflow']         = 'minimeta'
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if( params.notrim ){
    summary['Trimming Step'] = 'Skipped'
} else {
    summary["Trimming Step"] = 'Trim Glore'
}
if(workflow.profile == 'awsbatch'){
    summary['AWS Region']    = params.awsregion
    summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
    summary['E-mail Address']  = params.email
    summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
// 打印摘要
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// 生成工作流摘要 YAML 文件（用于 MultiQC）
def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-scgs-minimeta-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'gongyh/nf-core-scgs MINIMETA Workflow Summary'
    section_href: 'https://github.com/gongyh/nf-core-scgs'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v != null ? v : '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()
    return yaml_file
}

// ========================= 模块导入 ================================
// 导入 nf-core 标准模块
include { FASTQC                } from '../modules/nf-core/fastqc/main'   // 质量评估
include { MULTIQC               } from '../modules/nf-core/multiqc/main'  // 报告汇总

// 导入本地自定义模块
include { TRIMGALORE            } from '../modules/local/trimgalore'      // 质量修剪（Trim Galore）
include { CONCATENATE           } from '../modules/local/concatenate'     // 合并 FASTQ（透传）
include { SPADES                } from '../modules/local/spades'          // 子样本组装（已修改为输出校正 reads）
include { COLLECT_CORRECTED     } from '../modules/local/collect_corrected' // 收集校正 reads
include { SPADES_JOINT          } from '../modules/local/spades_joint'    // 联合组装（大内存）
include { RENAME_SUPERCONTIGS   } from '../modules/local/rename_supercontigs' // 重命名 super_contigs
include { BOWTIE2_INDEX         } from '../modules/local/bowtie2_index'   // 构建 Bowtie2 索引
include { BOWTIE2_ALIGN_SUBSAMPLE } from '../modules/local/bowtie2_align_subsample' // 子样本比对
include { SAMTOOLS_COVERAGE     } from '../modules/local/samtools_coverage' // 计算覆盖度
include { MERGE_COVERAGE        } from '../modules/local/merge_coverage'   // 合并覆盖度矩阵
include { COOCCURRENCE_BINNING  } from '../modules/local/cooccurrence_binning' // 共现分箱
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions/main' // 版本收集
include { OUTPUT_DOCUMENTATION  } from '../modules/local/output_documentation' // 文档输出

// 子工作流（用于邮件通知和完成摘要）
include { completionEmail       } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { completionSummary     } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'

def multiqc_report = []

workflow MINIMETA {
    // 输入通道：reads_ch，由上层传入或由 --reads 自动构建
    // 每个元素为 [meta, [r1, r2]]，meta.id 为子样本 ID
    take: reads_ch

    main:
    // 初始化版本收集通道
    ch_versions = Channel.empty()

    // 1. 合并 FASTQ（透传，统一命名）
    // 对应原始流程中的 concatenate 规则
    CONCATENATE( reads_ch )
    // 收集版本信息
    ch_versions = ch_versions.mix(CONCATENATE.out.versions)

    // 2. 原始数据质量评估（FastQC）
    FASTQC( CONCATENATE.out.merged_reads )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // 3. 质量修剪（Trim Galore）
    // 若用户指定 --notrim，则跳过修剪步骤
    if (params.notrim) {
        trimmed_reads = CONCATENATE.out.merged_reads
    } else {
        TRIMGALORE( CONCATENATE.out.merged_reads )
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
        trimmed_reads = TRIMGALORE.out.reads
    }

    // 4. 子样本组装（单细胞模式 SPAdes）
    // 该模块已修改，输出 contigs 和校正 reads (p1_corr, p2_corr, s_corr)
    SPADES( trimmed_reads )
    // 捕获输出
    contig      = SPADES.out.contig      // 每个子样本的 contigs（可选）
    p1_corr     = SPADES.out.p1_corr     // 校正后的 R1 reads
    p2_corr     = SPADES.out.p2_corr     // 校正后的 R2 reads
    s_corr      = SPADES.out.s_corr      // 校正后的单端 reads
    ch_versions = ch_versions.mix(SPADES.out.versions)

    // 5. 收集所有子样本的校正 reads 并合并 
    // 使用 collect() 将每个子样本的文件列表收集为单个列表
    p1_list = p1_corr.collect()
    p2_list = p2_corr.collect()
    s_list  = s_corr.collect()
    // 合并为三个总文件
    COLLECT_CORRECTED( p1_list, p2_list, s_list )
    ch_versions = ch_versions.mix(COLLECT_CORRECTED.out.versions)

    // 6. 联合组装（大内存 SPAdes）
    // 使用合并后的校正 reads 进行二次组装
    SPADES_JOINT( COLLECT_CORRECTED.out.r1, COLLECT_CORRECTED.out.r2, COLLECT_CORRECTED.out.s )
    ch_versions = ch_versions.mix(SPADES_JOINT.out.versions)

    // 7. 重命名 super_contigs 
    // 将 SPAdes 默认的 contig 名称（如 NODE_1_length_...）改为 SuperContig_1 等
    RENAME_SUPERCONTIGS( SPADES_JOINT.out.contigs )
    ch_versions = ch_versions.mix(RENAME_SUPERCONTIGS.out.versions)
    super_contigs = RENAME_SUPERCONTIGS.out.super_contigs

    //  8. 为 super_contigs 建立 Bowtie2 索引 
    BOWTIE2_INDEX( super_contigs )
    ch_versions = ch_versions.mix(BOWTIE2_INDEX.out.versions)
    index = BOWTIE2_INDEX.out.index

    //9. 每个子样本的 reads 比对到 super_contigs
    // 使用修剪后的 reads（trimmed_reads）进行比对
    BOWTIE2_ALIGN_SUBSAMPLE( trimmed_reads, index )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN_SUBSAMPLE.out.versions)
    bam_files = BOWTIE2_ALIGN_SUBSAMPLE.out.bam

    // 10. 计算覆盖度矩阵
    // 从 BAM 文件统计每个 contig 在每个子样本中的覆盖深度
    SAMTOOLS_COVERAGE( bam_files, super_contigs )
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions)
    cov_files = SAMTOOLS_COVERAGE.out.coverage

    // 11. 合并所有子样本的覆盖度为一个矩阵 
    // 输出为 TSV 文件，行是 contig，列是子样本，值为归一化覆盖度
    MERGE_COVERAGE( cov_files, super_contigs )
    ch_versions = ch_versions.mix(MERGE_COVERAGE.out.versions)
    coverage_matrix = MERGE_COVERAGE.out.matrix

    // 12. 共现分箱（核心算法）
    // 基于覆盖度矩阵进行：二值化 → Fisher 检验 → t‑SNE → DBSCAN 聚类
    COOCCURRENCE_BINNING( coverage_matrix, super_contigs )
    ch_versions = ch_versions.mix(COOCCURRENCE_BINNING.out.versions)
    clusters = COOCCURRENCE_BINNING.out.clusters   // contig 与簇对应表
    bins_dir = COOCCURRENCE_BINNING.out.bins       // 每个簇的 FASTA 序列

    // 13. 收集所有工具的版本信息
    GET_SOFTWARE_VERSIONS( ch_versions.unique().collectFile(name: 'collated_versions.yml') )
    ch_multiqc_versions = GET_SOFTWARE_VERSIONS.out.mqc_yml

    // 14. 准备 MultiQC 输入文件 
    workflow_summary = create_workflow_summary(summary)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    // 添加工作流摘要
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // 添加 FastQC 报告
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    // 添加 TrimGalore 日志
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]))
    // 添加版本信息
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_versions)

    // 运行 MultiQC 生成汇总报告
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

    // 生成输出文档（HTML）
    OUTPUT_DOCUMENTATION(ch_output_docs)

    // 15. 输出最终结果（供外部使用）
    emit:
    super_contigs = super_contigs
    coverage_matrix = coverage_matrix
    clusters = clusters
    bins = bins_dir
}

/*
 * Completion e-mail notification
 * 流程完成时的邮件通知
 */
workflow.onComplete {
    if (params.email){
        completionEmail(summary_params,
            params.email,
            null,
            false,
            params.outdir,
            log,
            multiqc_report.getVal()
        )
    }
    completionSummary()
}

// nf-core 风格的头部信息（彩色）
def nfcoreHeader(){
    // 定义颜色代码
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,'${c_reset}
    ${c_purple}  gongyh/nf-core-scgs MINIMETA v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}
