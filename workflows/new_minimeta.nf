def helpMessage() {
    log.info nfcoreHeader()
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
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// default values
params.single_end = false
params.notrim = false
params.saveTrimmed = false

custom_runName = workflow.runName
single_end = params.single_end

if(workflow.profile == 'awsbatch') {
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)
ch_multiqc_custom_config = Channel.empty()
ch_multiqc_logo = Channel.empty()
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

/*
 * Create a channel for input read files
 */
if(params.readPaths){
    if(single_end){
        read_files_fastqc = read_files_trimming =
        Channel.from(params.readPaths, checkIfExists: true)
            .map { row -> def meta=[:];
                    meta.id = row[0];
                    meta.single_end = single_end;
                    [meta, [file(row[1][0]), file(row[1][1])]]}
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
    } else {
        read_files_fastqc = read_files_trimming =
        Channel.from(params.readPaths)
            .map { row -> def meta=[:];
                    meta.id = row[0];
                    meta.single_end = single_end;
                    [meta, [file(row[1][0]), file(row[1][1])]]}
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
    }
} else {
    if (single_end) {
        read_files_fastqc = read_files_trimming =
        Channel.fromFilePairs(params.reads, size:1, checkIfExists: true)
            .map { it ->
                def meta = [:];
                meta.id = it[0].replaceFirst(~/\.[^\.]+$/, '');
                meta.single_end = single_end;
                [meta, [file(it[1][0])]]}

    } else {
        read_files_fastqc = read_files_trimming =
        Channel.fromFilePairs(params.reads, size:2, checkIfExists: true)
            .map { it ->
                def meta = [:];
                meta.id = it[0].replaceFirst(~/\.[^\.]+$/, '');
                meta.single_end = single_end;
                [meta, [file(it[1][0]), file(it[1][1])]]}
    }
}

// Header log info
log.info nfcoreHeader()
def summary = [:]
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
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

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

// Import modules
include { FASTQC                } from '../modules/nf-core/fastqc/main'
include { MULTIQC               } from '../modules/nf-core/multiqc/main'

include { TRIMGALORE            } from '../modules/local/trimgalore'
include { BBNORM                } from '../modules/local/bbnorm'
include { SPADES                } from '../modules/local/spades'
include { SPADES_JOINT          } from '../modules/local/spades_joint'
include { RENAME_SUPERCONTIGS   } from '../modules/local/rename_supercontigs'
include { BOWTIE2_REMAP         } from '../modules/local/bowtie2_remap'
include { BOWTIE2_ALIGN         } from '../modules/local/bowtie2_align'
include { SAMTOOLS              } from '../modules/local/samtools'
include { MERGE_COVERAGE        } from '../modules/local/merge_coverage'
include { COOCCURRENCE_BINNING  } from '../modules/local/cooccurrence_binning'
include { OUTPUT_DOCUMENTATION  } from '../modules/local/output_documentation'
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions/main'

/** subworkflow */
include { completionEmail       } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { completionSummary     } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'

// MULTIQC
def multiqc_report = []

workflow MINIMETA {
    ch_versions = Channel.empty()

    // FASTQC
    ch_multiqc_fastqc = Channel.empty()
    FASTQC ( read_files_fastqc )
    ch_versions       = ch_versions.mix(FASTQC.out.versions)
    ch_multiqc_fastqc = FASTQC.out.zip

    // TRIM_GALORE
    trimmed_reads = Channel.empty()
    ch_multiqc_trim_log = Channel.empty()
    ch_multiqc_trim_zip = Channel.empty()
    if (params.notrim) {
        trimmed_reads = read_files_trimming.map{name, reads -> reads}
    } else {
        TRIMGALORE ( read_files_trimming )
        ch_multiqc_trim_log = TRIMGALORE.out.log
        ch_multiqc_trim_zip = TRIMGALORE.out.zip
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
        trimmed_reads = TRIMGALORE.out.reads
    }

    // BBNORM
    BBNORM(trimmed_reads)
    normalized_reads = BBNORM.out.fastq
    ch_versions = ch_versions.mix(BBNORM.out.versions)

    // SPADES
    SPADES(normalized_reads)
    contig        = SPADES.out.contig
    contig_path   = SPADES.out.contig_path
    contig_graph  = SPADES.out.contig_graph
    ctg200        = SPADES.out.ctg200
    ctg           = SPADES.out.ctg
    assembly      = SPADES.out.assembly
    p1_corr       = SPADES_SS.out.p1_corr
    p2_corr       = SPADES_SS.out.p2_corr
    s_corr        = SPADES_SS.out.s_corr
    ch_versions   = ch_versions.mix(SPADES.out.versions)

    // ========== 第二阶段：手动联合组装 ==========
    // 收集所有子样本的校正 reads
    p1_ch = p1_corr.collect()
    p2_ch = p2_corr.collect()
    s_ch  = s_corr.collect()
    COLLECT_CORRECTED( p1_ch, p2_ch, s_ch )
    ch_versions = ch_versions.mix(COLLECT_CORRECTED.out.versions)

    // 2. 联合组装（大内存）
    p1_list = p1_corr.collect()
    p2_list = p2_corr.collect()
    s_list  = s_corr.collect()
    SPADES_JOINT( p1_list, p2_list, s_list )
    ch_versions = ch_versions.mix(SPADES_JOINT.out.versions)
    // 联合组装输出为 gzipped FASTA
    joint_contigs_gz = SPADES_JOINT.out.contigs

    // 解压联合组装结果
    process UNGZIP {
        input:
        path in_gz
        output:
        path "super_contigs.fasta"
        script:
        "gunzip -c ${in_gz} > super_contigs.fasta"
    }
    UNGZIP( joint_contigs_gz )
    super_contigs = UNGZIP.out

    // 3. 重命名 super_contigs
    RENAME_SUPERCONTIGS( SPADES_JOINT.out.contigs )
    ch_versions = ch_versions.mix(RENAME_SUPERCONTIGS.out.versions)
    super_contigs = RENAME_SUPERCONTIGS.out.super_contigs

    // 4. 建立 Bowtie2 索引
    // 注意：BOWTIE2_REMAP 输入要求 tuple val(meta), path(contigs)
    BOWTIE2_REMAP( [ [id:'ref'], final_super_contigs ] )
    ch_versions = ch_versions.mix(BOWTIE2_REMAP.out.versions)
    index_dir = BOWTIE2_REMAP.out.index

    // 5. 子样本 reads 比对到 super_contigs
    // trimmed_reads 格式为 [meta, [r1, r2]]；BOWTIE2_ALIGN 需要 (reads, index, save_unaligned, sort_bam)
    BOWTIE2_ALIGN( trimmed_reads, index_dir, false, true )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
    bam_files = BOWTIE2_ALIGN.out.bam

    // 6. 计算每个 contig 在每个子样本中的覆盖深度
    // 需修改 samtools.nf 增加 coverage.tsv 输出
    SAMTOOLS( bam_files, final_super_contigs )
    ch_versions = ch_versions.mix(SAMTOOLS.out.versions)
    cov_files = SAMTOOLS.out.coverage

    // 7. 合并所有子样本的覆盖度矩阵（行=contig，列=子样本，值=归一化覆盖度）
    MERGE_COVERAGE( cov_files.collect(), final_super_contigs )
    ch_versions = ch_versions.mix(MERGE_COVERAGE.out.versions)
    coverage_matrix = MERGE_COVERAGE.out.matrix

    // 8. 共现分箱（Fisher检验 + t‑SNE + DBSCAN）
    COOCCURRENCE_BINNING( coverage_matrix, final_super_contigs )
    ch_versions = ch_versions.mix(COOCCURRENCE_BINNING.out.versions)
    clusters = COOCCURRENCE_BINNING.out.clusters
    bins_dir = COOCCURRENCE_BINNING.out.bins

    // GET_SOFTWARE_VERSIONS
    ch_multiqc_versions = Channel.empty()
    GET_SOFTWARE_VERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    ch_multiqc_versions = GET_SOFTWARE_VERSIONS.out.mqc_yml

    // MODULE: MULTIQC
    workflow_summary = create_workflow_summary(summary)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_fastqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_trim_log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_trim_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_versions)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

    OUTPUT_DOCUMENTATION(ch_output_docs)
}

/*
 * Completion e-mail notification
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

def nfcoreHeader(){
    // Log colors ANSI codes
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
