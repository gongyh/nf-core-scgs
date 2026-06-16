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

    Database options:
    --mmseqs_db                   Path to the mmseqs database

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
// default values
params.single_end = false
params.notrim = false
params.saveTrimmed = false
params.mmseqs_db = null
custom_runName = workflow.runName
single_end = params.single_end

if(workflow.profile == 'awsbatch') {
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Configure Checkm2 database
checkm2_db = false
if (params.checkm2_db) {
    checkm2_db  = file(params.checkm2_db)
    if ( !checkm2_db.exists() ) exit 1, "CheckM2 database not found: ${params.checkm2_db}"
} else {
    checkm2_db = file("/dev/null")
}

//kofam database
kofam_profile = false
if (params.kofam_profile) {
    kofam_profile = file(params.kofam_profile)
    if( !kofam_profile.exists() ) exit 1, "KOfam profile database not found: ${params.kofam_profile}"
} else {
    kofam_profile = file("/dev/null")
}

kofam_kolist = false
if (params.kofam_kolist) {
    kofam_kolist = file(params.kofam_kolist)
    if( !kofam_kolist.exists() ) exit 1, "KOfam ko_list file not found: ${params.kofam_kolist}"
} else {
    kofam_kolist = file("/dev/null")
}

//eggnog database
eggnog_db = false
if (params.eggnog_db) {
    eggnog_db = file(params.eggnog_db)
    if( !eggnog_db.exists() ) exit 1, "EggNOG database not found: ${params.eggnog_db}"
} else {
    eggnog_db = file("/dev/null")
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
include { SPADES as SPADES_JOINT            } from '../modules/local/spades'
include { BOWTIE2_REMAP                     } from '../modules/local/bowtie2_remap'
include { REMAP                             } from '../modules/local/remap'
include { MERGE_BAMS                        } from '../modules/local/merge_bams'
include { SAMTOOLS_FAIDX                    } from '../modules/local/samtools_faidx'
include { PREPARE_FEATURES_SINGLE           } from '../subworkflows/local/prepare_features_single'
include { PREPARE_FEATURES_MULTI            } from '../subworkflows/local/prepare_features_multi'
include { COOCCURRENCE_BINNING              } from '../modules/local/binning'
include { EXTRACT_BINS                      } from '../modules/local/extract_bins'
include { MMSEQS_CONTIG_TAXONOMY            } from '../subworkflows/local/mmseqs_contig_taxonomy'
include { SEMIBIN2                          } from '../modules/local/semibin2'
include { TAXVAMB_INTEGRATION               } from '../subworkflows/local/taxvamb_integration'
include { DAS_TOOL                          } from '../modules/local/das_tool'
include { CHECKM2                           } from '../modules/local/checkm2'
include { PROKKA                            } from '../modules/local/prokka'
include { KOFAMSCAN                         } from '../modules/local/kofamscan'
include { EGGNOG                            } from '../modules/local/eggnog'
include { OUTPUT_DOCUMENTATION              } from '../modules/local/output_documentation'
include { GET_SOFTWARE_VERSIONS             } from '../modules/local/get_software_versions/main'

// MULTIQC
def multiqc_report = []

workflow MINIMETA {
    main:
    display_header()
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
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
    // Sort
    p1_list = corrected_reads.map { meta, reads -> reads[0] }.collect()
    p2_list = corrected_reads.map { meta, reads -> reads[1] }.collect()

    //Merge_corrected
    MERGE_CORRECTED( p1_list, p2_list )
    joint_reads = MERGE_CORRECTED.out.r1
        .combine(MERGE_CORRECTED.out.r2)
        .map { r1, r2 -> [ [id:'merged', single_end:false], [r1, r2] ] }

    //SPADES_JOINT
    SPADES_JOINT( joint_reads )
    ch_versions = ch_versions.mix(SPADES_JOINT.out.versions)

    //BOWTIE2_REMAP
    BOWTIE2_REMAP( SPADES_JOINT.out.contig )
    ch_versions = ch_versions.mix(BOWTIE2_REMAP.out.versions)
    //REMAP
    remap_input = trimmed_reads.combine(BOWTIE2_REMAP.out.index).map {
        [it[0] + [id_index: 'merged'], it[1], it[3]]
    }
    REMAP(remap_input, params.allow_multi_align)
    ch_versions = ch_versions.mix(REMAP.out.versions)

    //MERGE BAMS
    ch_bam_list = REMAP.out.bam.map{ meta, bam -> bam }.collect()
    MERGE_BAMS( ch_bam_list )
    ch_merged_bam = MERGE_BAMS.out.merged_bam
    ch_bam_for_coverage = ch_merged_bam.map { bam -> [ [id:'merged'], bam, [] ] }

    //PREPARE_FEATURES
    ch_fasta = SPADES_JOINT.out.contig
    SAMTOOLS_FAIDX( ch_fasta )
    ch_fai = SAMTOOLS_FAIDX.out.fai
    PREPARE_FEATURES_SINGLE( ch_fasta, ch_fai, ch_bam_for_coverage )
    ch_single_coverage = PREPARE_FEATURES_SINGLE.out.coverage_matrix
    PREPARE_FEATURES_MULTI( ch_fasta, ch_fai, REMAP.out.bam )
    ch_multi_coverage = PREPARE_FEATURES_MULTI.out.coverage_matrix

    // binning
    ch_assembly = SPADES_JOINT.out.contig.map { it[1] }
    ch_all_s2b = Channel.empty()
    ch_versions = Channel.empty()

    //COOCCURRENCE
    COOCCURRENCE_BINNING( ch_multi_coverage )
    EXTRACT_BINS(COOCCURRENCE_BINNING.out.clusters, ch_assembly)
    ch_all_s2b = ch_all_s2b.mix( EXTRACT_BINS.out.scaffolds2bin.map { file -> ['COOCCURRENCE', file] } )
    ch_versions = ch_versions.mix( COOCCURRENCE_BINNING.out.versions )

    //SEMIBIN2
    SEMIBIN2(ch_assembly, ch_merged_bam)
    ch_semibin2_s2b = SEMIBIN2.out.scaffolds2bin
        .map { file -> ['SEMIBIN2', file] }
        .filter { it[1].size() > 0 }
    ch_all_s2b = ch_all_s2b.mix(ch_semibin2_s2b)
    ch_versions = ch_versions.mix( SEMIBIN2.out.versions )

    if (params.mmseqs_db ) {
        //MMseqs_TAXA
        ch_mmseqs_db = channel.fromPath( params.mmseqs_db )
        MMSEQS_CONTIG_TAXONOMY( ch_assembly, ch_mmseqs_db )
        //SEMIBIN2_Semi
    }

    // TaxVAMB
    TAXVAMB_INTEGRATION( ch_assembly, ch_single_coverage )
    ch_all_s2b = ch_all_s2b.mix( TAXVAMB_INTEGRATION.out.scaffolds2bin.map { file -> ['TAXVAMB', file] } )
    ch_versions = ch_versions.mix( TAXVAMB_INTEGRATION.out.versions )

    // DAS TOOL
    ch_s2b_list = ch_all_s2b.flatten().toList()
    DAS_TOOL(ch_assembly, ch_s2b_list)
    ch_bins_dir = DAS_TOOL.out.bins
    ch_versions = ch_versions.mix(DAS_TOOL.out.versions)

    // CHECKM2
    CHECKM2(ch_bins_dir, "fa", file(params.checkm2_db ?: "/dev/null"))
    ch_versions = ch_versions.mix(CHECKM2.out.versions)
    ch_multiqc_checkm2 = CHECKM2.out.mqc_tsv
    ch_multiqc_files = ch_multiqc_files.mix(CHECKM2.out.mqc_tsv)
    //
    ch_bins_for_prokka = ch_bins_dir.flatMap { bin_dir ->
        def bin_files = file(bin_dir).listFiles().findAll { it.name.endsWith('.fa') }
        if (!bin_files) {
            log.warn "No .fa files found in ${bin_dir}, skipping PROKKA"
            return []
        }
        bin_files.collect { bin_file ->
            [ [id: bin_file.baseName], bin_file ]
        }
    }

    //PROKKA
    PROKKA(ch_bins_for_prokka, [])
    ch_versions = ch_versions.mix(PROKKA.out.versions)

    // KOFAMSCAN
    if (params.kofam) {
        KOFAMSCAN(PROKKA.out.faa, kofam_profile, kofam_kolist)
        ch_versions = ch_versions.mix(KOFAMSCAN.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(KOFAMSCAN.out.kofamscan.collect().ifEmpty([]))
    }

    // EGGNOG
    if (params.eggnog) {
        EGGNOG(PROKKA.out.faa, eggnog_db)
        ch_versions = ch_versions.mix(EGGNOG.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(EGGNOG.out.annotations.collect().ifEmpty([]))
    }

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
    ch_multiqc_files = ch_multiqc_files.mix(SPADES_JOINT.out.mqc_tsv.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(REMAP.out.mqc_tsv.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_FEATURES_MULTI.out.coverage_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_FEATURES_SINGLE.out.coverage_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(COOCCURRENCE_BINNING.out.mqc_tsv.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(EXTRACT_BINS.out.mqc_tsv.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SEMIBIN2.out.mqc_tsv.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix( TAXVAMB_INTEGRATION.out.mqc_tsv.ifEmpty([]) )
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

