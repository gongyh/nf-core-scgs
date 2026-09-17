def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    nextflow run gongyh/nf-core-scgs --reads '*_R{1,2}.fastq.gz' --minimeta -profile docker

    Workflow selection:
    --minimeta                    Run the MINIMETA workflow

    Input options:
    --reads <glob>                Input reads glob (default: data/*{1,2}.fastq.gz)
    --readPaths <list>            Structured sample/read list supplied in a Nextflow config
    --single_end                  Treat input reads as single-end

    Read processing:
    --notrim                      Skip adapter and quality trimming
    --saveTrimmed                 Publish trimmed reads
    --clip_r1 <int>               Remove bases from the 5' end of read 1
    --clip_r2 <int>               Remove bases from the 5' end of read 2
    --three_prime_clip_r1 <int>   Remove bases from the 3' end of read 1 after trimming
    --three_prime_clip_r2 <int>   Remove bases from the 3' end of read 2 after trimming
    --allow_multi_align           Retain secondary and unmapped remapping alignments

    Binning and quality assessment:
    --min_length <int>            Minimum contig length for co-occurrence binning (default: 10000)
    --cooccurrence_eps <number>   Distance threshold for co-occurrence binning (default: 0.05)
    --run_cooccurrence_checkm     Run CheckM2 on co-occurrence bins when --checkm2_db is available
    --checkm2_db <path>           CheckM2 database
    --mmseqs_db <path>            MMseqs2 database for contig taxonomy and SemiBin2
    --metabuli_db <path>          MetaBuli database; enables TaxVAMB integration
    --DNABERTS_dir <path>         DNABERT-S model directory; enables DCVBIN integration

    Functional annotation:
    --kofam                       Run KOfam annotation when profile and KO-list files are available
    --kofam_profile <path>        KOfam profile database
    --kofam_kolist <path>         KOfam KO-list file
    --eggnog                      Run EggNOG annotation when --eggnog_db is available
    --eggnog_db <path>            EggNOG database

    Output and execution:
    --outdir <path>               Output directory (default: ./results)
    --multiqc_config <path>       Custom MultiQC configuration file
    --email <address>             Address for the completion email
    --maxMultiqcEmailFileSize     Maximum MultiQC email attachment size in bytes (default: 25 MB)
    --monochrome_logs             Disable coloured log output
    --help                        Display this help message
    --awsqueue <name>             AWS Batch job queue
    --awsregion <region>          AWS Batch region
    -profile                      Configuration profile(s), for example: docker, singularity, conda
    """.stripIndent()
}

def display_header(summary, custom_runName, single_end) {
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
include { FILTER_ASSEMBLY                   } from '../modules/local/filter_assembly'
include { COOCCURRENCE_BINNING              } from '../modules/local/binning'
include { CHECKM2 as CHECKM2_COOCCURRENCE   } from '../modules/local/checkm2'
include { EXTRACT_BINS                      } from '../modules/local/extract_bins'
include { SEMIBIN2                          } from '../modules/local/semibin2'
include { MMSEQS_CONTIG_TAXONOMY            } from '../subworkflows/local/mmseqs_contig_taxonomy'
include { MMSEQS2SEMIBIN                    } from '../modules/local/mmseqs2semibin'
include { TAXVAMB_INTEGRATION               } from '../subworkflows/local/taxvamb_integration'
include { FILTER_CONTIGS                    } from '../modules/local/filter_contigs'
include { FILTER_BAM                        } from '../modules/local/filter_bam'
include { DCVBIN                            } from '../subworkflows/local/dcvbin'
include { DAS_TOOL                          } from '../modules/local/das_tool'
include { CHECKM2                           } from '../modules/local/checkm2'
include { PROKKA                            } from '../modules/local/prokka'
include { KOFAMSCAN                         } from '../modules/local/kofamscan'
include { EGGNOG                            } from '../modules/local/eggnog'
include { OUTPUT_DOCUMENTATION              } from '../modules/local/output_documentation'
include { GET_SOFTWARE_VERSIONS             } from '../modules/local/get_software_versions/main'

workflow MINIMETA {
    main:
    /*
 * SET UP CONFIGURATION VARIABLES
 */
// default values
params.reads = "data/*{1,2}.fastq.gz"
params.outdir = "./results"
params.notrim = false
params.awsregion = "eu-west-1"
params.awsqueue = "default"
params.config_profile_description = null
params.config_profile_contact = null
params.config_profile_url = null
params.email = null
params.maxMultiqcEmailFileSize = 25 * 1024 * 1024
params.single_end = false
params.checkm2_db = null
params.kofam_profile = null
params.kofam_kolist = null
params.eggnog_db = null
params.multiqc_config = "$baseDir/assets/multiqc_config.yml"
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0
params.readPaths = null
params.saveTrimmed = false
params.bulk = false
params.mg = false
params.allow_multi_align = false
params.min_length = 10000
params.run_cooccurrence_checkm = false
params.cooccurrence_eps = 0.05
params.mmseqs_db = null
params.metabuli_db = null
params.DNABERTS_dir = null
params.kofam = true
params.eggnog = true
params.monochrome_logs = false
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
ch_multiqc_config = channel.fromPath(params.multiqc_config, checkIfExists: true)
ch_multiqc_custom_config = channel.empty()
ch_multiqc_logo = channel.empty()
ch_output_docs = channel.fromPath("$baseDir/docs/output.md")


/*
 * Create a channel for input read files
 */
if(params.readPaths){
    if(single_end){
        read_files_fastqc = channel.from(params.readPaths, checkIfExists: false)
            .map { row -> def meta=[:];
                    meta.id = row[0];
                    meta.single_end = single_end;
                    [meta, [file(row[1][0]), file(row[1][1])]]}
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        read_files_trimming = read_files_fastqc
    } else {
        read_files_fastqc = channel.from(params.readPaths)
            .map { row -> def meta=[:];
                    meta.id = row[0];
                    meta.single_end = single_end;
                    [meta, [file(row[1][0]), file(row[1][1])]]}
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
        read_files_trimming = read_files_fastqc
    }
} else {
    if (single_end) {
        read_files_fastqc = channel.fromFilePairs(params.reads, size:1, checkIfExists: false)
            .map { it ->
                def meta = [:];
                meta.id = it[0].replaceFirst(~/\.[^\.]+$/, '');
                meta.single_end = single_end;
                [meta, [file(it[1][0])]]}
        read_files_trimming = read_files_fastqc

    } else {
        read_files_fastqc = channel.fromFilePairs(params.reads, size:2, checkIfExists: false)
            .map { it ->
                def meta = [:];
                meta.id = it[0].replaceFirst(~/\.[^\.]+$/, '');
                meta.single_end = single_end;
                [meta, [file(it[1][0]), file(it[1][1])]]}
        read_files_trimming = read_files_fastqc
    }
}

summary = [:]


    display_header(summary, custom_runName, single_end)
    ch_published = channel.empty()
    ch_multiqc_files = channel.empty()
    // FASTQC
    ch_multiqc_fastqc = channel.empty()
    FASTQC ( read_files_fastqc )
    ch_vendor_versions = FASTQC.out.versions
    ch_multiqc_fastqc = FASTQC.out.zip
    ch_published = ch_published.mix(FASTQC.out.html.map { result -> [destination: 'fastqc', files: result] })
    ch_published = ch_published.mix(FASTQC.out.zip.map { result -> [destination: 'fastqc/zips', files: result] })

    // TRIM_GALORE
    trimmed_reads = channel.empty()
    ch_multiqc_trim_log = channel.empty()
    ch_multiqc_trim_zip = channel.empty()
    if (params.notrim) {
        trimmed_reads = read_files_trimming
    } else {
        trimgalore = TRIMGALORE(read_files_trimming)
        ch_multiqc_trim_log = trimgalore.map { result -> tuple(result.meta, result.logs) }
        ch_multiqc_trim_zip = trimgalore.map { result -> tuple(result.meta, result.fastqc) }
        trimmed_reads = trimgalore.map { result ->
            def reads = result.meta.single_end ? [result.single_read] : [result.read1, result.read2]
            tuple(result.meta, reads)
        }
        ch_published = ch_published.mix(trimgalore.map { result -> [destination: "trim_galore/${result.meta.id}", files: result.fastqc] })
        ch_published = ch_published.mix(trimgalore.map { result -> [destination: "trim_galore/${result.meta.id}", files: result.logs] })
        if (params.saveTrimmed) {
            ch_published = ch_published.mix(trimmed_reads.map { _meta, reads -> [destination: 'trim_galore', files: reads] })
        }
    }

    // BBNORM
    bbnorm = BBNORM(trimmed_reads)
    normalized_reads = bbnorm.map { result ->
        def reads = result.meta.single_end ? [result.single_fastq] : [result.fastq1, result.fastq2]
        tuple(result.meta, reads)
    }

    // Performs read error correction for each minimeta sample
    read_correction = READ_CORRECTION(normalized_reads.map { meta, reads ->
        def meta_clone = meta.clone()
        meta_clone.only_error_correction = true;
        tuple(meta_clone, reads)
    })
    corrected_reads = read_correction.map { result ->
        def reads = result.meta.single_end ? [result.corrected_read] : [result.corrected_read, result.corrected_read2]
        tuple(result.meta, reads)
    }
    ch_published = ch_published.mix(read_correction.map { result -> [destination: 'spades', files: result] })
    // Sort
    p1_list = corrected_reads.map { meta, reads -> reads[0] }.collect()
    p2_list = corrected_reads.map { meta, reads -> reads[1] }.collect()

    //Merge_corrected
    merge_corrected = MERGE_CORRECTED(p1_list, p2_list)
    joint_reads = merge_corrected.map { result -> [ [id:'merged', single_end:false], [result.r1, result.r2] ] }
    ch_published = ch_published.mix(merge_corrected.map { result -> [destination: 'merged', files: result] })

    //SPADES_JOINT
    spades_joint = SPADES_JOINT(joint_reads)
    ch_published = ch_published.mix(spades_joint.map { result -> [destination: 'spades', files: result] })

    //BOWTIE2_REMAP
    ch_spades_contig = spades_joint.map { result -> tuple(result.meta, result.contig) }
    bowtie2_remap = BOWTIE2_REMAP(ch_spades_contig)
    //REMAP
    remap_input = trimmed_reads.combine(bowtie2_remap.map { result -> tuple(result.meta, result.index) }).map { entry ->
        tuple(entry[0] + [id_index: 'merged'], entry[1], entry[3])
    }
    remap = REMAP(remap_input, params.allow_multi_align)
    ch_published = ch_published.mix(remap.map { result -> [destination: 'remap', files: result] })

    //MERGE BAMS
    ch_remap_bam = remap.map { result -> tuple(result.meta, result.bam) }
    ch_bam_list = ch_remap_bam.map { _meta, bam -> bam }.collect()
    merge_bams = MERGE_BAMS(ch_bam_list)
    ch_merged_bam = merge_bams.map { result -> result.merged_bam }
    ch_bam_for_coverage = ch_merged_bam.map { bam -> tuple([id: 'merged'], bam, [] as List<Path>) }
    ch_published = ch_published.mix(merge_bams.map { result -> [destination: 'merged_bam', files: result] })

    //PREPARE_FEATURES
    ch_fasta = ch_spades_contig
    samtools_faidx = SAMTOOLS_FAIDX(ch_fasta)
    ch_fai = samtools_faidx.map { result -> [result.meta, result.fai] }
    PREPARE_FEATURES_SINGLE( ch_fasta, ch_fai, ch_bam_for_coverage )
    ch_single_coverage = PREPARE_FEATURES_SINGLE.out.coverage_matrix
    ch_published = ch_published.mix(PREPARE_FEATURES_SINGLE.out.published)
    PREPARE_FEATURES_MULTI( ch_fasta, ch_fai, ch_remap_bam )
    ch_multi_coverage = PREPARE_FEATURES_MULTI.out.coverage_matrix
    ch_published = ch_published.mix(PREPARE_FEATURES_MULTI.out.published)

    // binning
    ch_assembly = spades_joint.map { result -> result.contig }
    ch_all_s2b = channel.empty()

    def min_len = params.min_length ?: 10000
    filter_assembly = FILTER_ASSEMBLY(ch_assembly, min_len)
    ch_filtered_assembly = filter_assembly.map { result -> result.filtered }
    //COOCCURRENCE
    ch_filtered_ids = filter_assembly.map { result -> result.filtered_ids }
    cooccurrence_binning = COOCCURRENCE_BINNING(ch_multi_coverage, ch_filtered_ids)
    extract_bins = EXTRACT_BINS(cooccurrence_binning.map { result -> result.clusters }, ch_assembly)
    ch_all_s2b = ch_all_s2b.mix(extract_bins.map { result -> ['COOCCURRENCE', result.scaffolds2bin] })
    ch_published = ch_published.mix(cooccurrence_binning.map { result -> [destination: 'cooccurrence_bins', files: result] })
    ch_published = ch_published.mix(extract_bins.map { result -> [destination: 'extracted_bins', files: result] })

    //
    if ( params.run_cooccurrence_checkm ) {
        if ( params.checkm2_db ) {
            checkm2_cooccurrence = CHECKM2_COOCCURRENCE(extract_bins.map { result -> result.bins }, 'fa', file(params.checkm2_db))
            ch_multiqc_files = ch_multiqc_files.mix(checkm2_cooccurrence.map { result -> result.mqc_tsv }.collect().ifEmpty([]))
            ch_published = ch_published.mix(checkm2_cooccurrence.map { result -> [destination: 'CheckM2', files: result] })
        } else {
            log.info "INFO: --run_cooccurrence_checkm is set, but --checkm2_db is not provided. Skipping CheckM2 for COOCCURRENCE."
        }
    }
    //MMseqs2
    if (params.mmseqs_db ) {
        //MMseqs_TAXA
        ch_mmseqs_input = ch_assembly.map { fasta -> [ [id: fasta.baseName], fasta ] }
        ch_mmseqs_db = channel.fromPath( params.mmseqs_db )

        MMSEQS_CONTIG_TAXONOMY( ch_mmseqs_input, ch_mmseqs_db )
        ch_mmseqs_taxonomy = MMSEQS_CONTIG_TAXONOMY.out.taxonomy
        ch_multiqc_files = ch_multiqc_files.mix(ch_mmseqs_taxonomy.collect().ifEmpty([]))
        ch_published = ch_published.mix(MMSEQS_CONTIG_TAXONOMY.out.published)

        mmseqs2semibin = MMSEQS2SEMIBIN(ch_mmseqs_taxonomy)
        ch_semibin_tax = mmseqs2semibin.map { result -> result.tax }
        //SEMIBIN2_Semi
        semibin2 = SEMIBIN2(ch_assembly, ch_merged_bam, ch_semibin_tax)
    } else {
        semibin2 = SEMIBIN2(ch_assembly, ch_merged_bam, null)
    }

    //SEMIBIN2
    ch_semibin2_s2b = semibin2
        .map { result -> ['SEMIBIN2', result.scaffolds2bin] }
        .filter { entry -> entry[1].size() > 0 }
    ch_all_s2b = ch_all_s2b.mix(ch_semibin2_s2b)
    ch_published = ch_published.mix(semibin2.map { result -> [destination: 'semibin2_bins', files: result] })

    if (params.metabuli_db) {
        taxvamb = TAXVAMB_INTEGRATION(ch_assembly, ch_single_coverage)
        ch_all_s2b = ch_all_s2b.mix(taxvamb.scaffolds2bin.map { _meta, file -> ['TAXVAMB', file] })
        ch_published = ch_published.mix(taxvamb.published)
        ch_taxvamb_mqc = taxvamb.mqc_tsv
    } else {
        ch_taxvamb_mqc = channel.empty()
    }

    if (params.DNABERTS_dir != null){
        //FILTERED
        ch_merged_bai = ch_merged_bam.map { bam -> file("${bam}.bai") }
        filter_contigs = FILTER_CONTIGS(ch_assembly, 2000)
        ch_filtered_fasta_with_meta = filter_contigs.map { result ->
            def fasta = result.filtered
            [ [id: fasta.baseName], fasta ]
        }
        filter_bam = FILTER_BAM(ch_filtered_fasta_with_meta, ch_merged_bam, ch_merged_bai)
        ch_filtered_bam = filter_bam.map { result -> result.filtered_bam }
        ch_published = ch_published.mix(filter_contigs.map { result -> [destination: 'filtered_fasta', files: result] })
        ch_published = ch_published.mix(filter_bam.map { result -> [destination: 'filtered_bam', files: result] })
        ch_bam_path = ch_filtered_bam
        //DCVBIN
        dcvbin = DCVBIN(ch_filtered_fasta_with_meta, ch_bam_path)
        ch_all_s2b = ch_all_s2b.mix(dcvbin.scaffolds2bin.map { _meta, file -> ['DCVBIN', file] })
        ch_multiqc_files = ch_multiqc_files.mix(dcvbin.mqc_tsv.ifEmpty([]))
        ch_published = ch_published.mix(dcvbin.published)
    }

    // DAS TOOL
    ch_s2b_list = ch_all_s2b.flatten().toList()
    das_tool = DAS_TOOL(ch_assembly, ch_s2b_list)
    ch_bins_dir = das_tool.map { result -> result.bins }
    ch_published = ch_published.mix(das_tool.map { result -> [destination: 'binning/das_tool', files: result] })

    // CHECKM2
    if (params.checkm2_db) {
        checkm2 = CHECKM2(ch_bins_dir, 'fa', checkm2_db)
        ch_multiqc_checkm2 = checkm2.map { result -> result.mqc_tsv }
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_checkm2)
        ch_published = ch_published.mix(checkm2.map { result -> [destination: 'CheckM2', files: result] })
    } else {
        ch_multiqc_checkm2 = channel.empty()
    }
    //
    ch_bins_for_prokka = ch_bins_dir.flatMap { bin_dir ->
        def bin_files = file(bin_dir).listFiles().findAll { entry -> entry.name.endsWith('.fa') }
        if (!bin_files) {
            log.warn "No .fa files found in ${bin_dir}, skipping PROKKA"
            return []
        }
        bin_files.collect { bin_file ->
            [ [id: bin_file.baseName], bin_file ]
        }
    }

    //PROKKA
    prokka = PROKKA(ch_bins_for_prokka, [] as List<Path>)
    ch_published = ch_published.mix(prokka.map { result -> [destination: 'prokka', files: result] })

    // KOFAMSCAN
    if (params.kofam && params.kofam_profile && params.kofam_kolist) {
        kofamscan = KOFAMSCAN(prokka.map { result -> tuple(result.meta, result.faa) }, kofam_profile, kofam_kolist)
        ch_multiqc_files = ch_multiqc_files.mix(kofamscan.map { result -> result.kofamscan }.collect().ifEmpty([]))
        ch_published = ch_published.mix(kofamscan.map { result -> [destination: 'kofam', files: result] })
    }

    // EGGNOG
    if (params.eggnog && params.eggnog_db) {
        eggnog = EGGNOG(prokka.map { result -> tuple(result.meta, result.faa) }, eggnog_db)
        ch_multiqc_files = ch_multiqc_files.mix(eggnog.map { result -> result.annotations }.collect().ifEmpty([]))
        ch_published = ch_published.mix(eggnog.map { result -> [destination: 'eggnog', files: result] })
    }

    // GET_SOFTWARE_VERSIONS
    ch_multiqc_versions = channel.empty()
    ch_local_versions = trimgalore.map { result -> result.versions }
        .mix(bbnorm.map { result -> result.versions })
        .mix(read_correction.map { result -> result.versions })
        .mix(spades_joint.map { result -> result.versions })
        .mix(bowtie2_remap.map { result -> result.versions })
        .mix(remap.map { result -> result.versions })
        .mix(merge_bams.map { result -> result.versions })
        .mix(samtools_faidx.map { result -> result.versions })
        .mix(PREPARE_FEATURES_SINGLE.out.versions)
        .mix(PREPARE_FEATURES_MULTI.out.versions)
        .mix(filter_assembly.map { result -> result.versions })
        .mix(cooccurrence_binning.map { result -> result.versions })
        .mix(extract_bins.map { result -> result.versions })
        .mix(semibin2.map { result -> result.versions })
        .mix(das_tool.map { result -> result.versions })
        .mix(prokka.map { result -> result.versions })
    software_versions = GET_SOFTWARE_VERSIONS(
        ch_local_versions
            .mix(ch_vendor_versions)
            .map { version ->
                def lines = version.text.readLines()
                def first_content = lines.find { line -> line.trim() && line.trim() != 'END_VERSIONS' }
                def indentation = first_content ? first_content.length() - first_content.stripLeading().length() : 0
                lines
                    .findAll { line -> line.trim() != 'END_VERSIONS' }
                    .collect { line -> indentation > 0 && line.length() >= indentation ? line.substring(indentation) : line }
                    .join('\n') + '\n'
            }
            .unique()
            .collectFile(name: 'collated_versions.yml', newLine: true)
    )
    ch_multiqc_versions = software_versions.map { result -> result.mqc_yml }
    ch_published = ch_published.mix(software_versions.map { result -> [destination: 'pipeline_info', files: [result.yml, result.mqc_yml]] })

    // MODULE: MULTIQC
    workflow_summary = create_workflow_summary(summary)
    ch_workflow_summary = channel.value(workflow_summary)

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_fastqc.collect { entry -> entry[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_trim_log.collect { entry -> entry[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_trim_zip.collect { entry -> entry[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(spades_joint.map { result -> result.mqc_tsv }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(remap.map { result -> result.mqc_tsv }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_FEATURES_MULTI.out.coverage_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_FEATURES_SINGLE.out.coverage_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(cooccurrence_binning.map { result -> result.mqc_tsv }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(extract_bins.map { result -> result.mqc_tsv }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(semibin2.map { result -> result.mqc_tsv }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_taxvamb_mqc.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_versions)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_published = ch_published.mix(MULTIQC.out.report.map { result -> [destination: 'MultiQC', files: result] })
    ch_published = ch_published.mix(MULTIQC.out.data.map { result -> [destination: 'MultiQC', files: result] })
    ch_published = ch_published.mix(MULTIQC.out.plots.map { result -> [destination: 'MultiQC', files: result] })
    ch_published = ch_published.mix(MULTIQC.out.versions.map { result -> [destination: 'MultiQC', files: result] })
    OUTPUT_DOCUMENTATION(ch_output_docs)

    emit:
    summary_params = channel.value(summary)
    multiqc_report = MULTIQC.out.report.toList()
    published = ch_published
}

def nfcoreHeader(){
    // Log colors ANSI codes
    def c_reset = params.monochrome_logs ? '' : "\033[0m";
    def c_dim = params.monochrome_logs ? '' : "\033[2m";
    def c_black = params.monochrome_logs ? '' : "\033[0;30m";
    def c_green = params.monochrome_logs ? '' : "\033[0;32m";
    def c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    def c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    def c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    def c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    def c_white = params.monochrome_logs ? '' : "\033[0;37m";

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
