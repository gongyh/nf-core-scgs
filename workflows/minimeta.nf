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

summary = [:]

def display_header() {
    // Header log info
    log.info nfcoreHeader()
    //def summary = [:]
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
}

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
include { FASTQC                            } from '../modules/nf-core/fastqc/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'

include { TRIMGALORE                        } from '../modules/local/trimgalore'
include { BBNORM                            } from '../modules/local/bbnorm'
include { SPADES as READ_CORRECTION; SPADES } from '../modules/local/spades'
include { MERGE_CORRECTED                   } from '../modules/local/merge_corrected'
include { OUTPUT_DOCUMENTATION              } from '../modules/local/output_documentation'
include { GET_SOFTWARE_VERSIONS             } from '../modules/local/get_software_versions/main'

// MULTIQC
def multiqc_report = []

workflow MINIMETA {
    main:
    display_header()
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

    // Performs read error correction for each minimeta sample
    READ_CORRECTION(normalized_reads.map { meta, reads ->
        def meta_clone = meta.clone()
        meta_clone.only_error_correction = true;
        [meta_clone, reads]
    })
    corrected_reads = READ_CORRECTION.out.reads
    ch_versions = ch_versions.mix(READ_CORRECTION.out.versions)
    //Merge_corrected
    p1_list = corrected_reads.map { meta, reads -> reads[0] }.collect()
    p2_list = corrected_reads.map { meta, reads -> reads[1] }.collect()

    MERGE_CORRECTED( p1_list, p2_list )
    ch_versions = ch_versions.mix(MERGE_CORRECTED.out.versions)

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
