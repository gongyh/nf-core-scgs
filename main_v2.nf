#!/usr/bin/env nextflow
/*
========================================================================================
                        gongyh/nf-core-scgs
========================================================================================
    gongyh/nf-core-scgs Analysis Pipeline.
    #### Homepage / Documentation
    https://github.com/gongyh/nf-core-scgs
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl=2

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run gongyh/nf-core-scgs --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
    --reads                       Path to input data (must be surrounded with quotes)
    -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    Options:
    --vcf                         Variantion graph construction
    --bulk                        WGS of bulk DNA, not MDA
    --nanopore                    Nanopore sequencing data, force single_end, assembly using Canu
    --genome                      Name of iGenomes reference
    --single_end                  Specifies that the input is single end reads
    --snv                         Enable detection of single nucleotide variation
    --cnv                         Enable detection of copy number variation
    --remap                       Remap trimmed reads to contigs
    --saturation                  Enable sequencing saturation analysis
    --ass                         Assemble using SPAdes
    --split                       Split the draft genomes and annotation
    --split_bac_level             Level of split for Bacteria
    --split_euk_level             Level of split for Eukaryota

    References:                     If not specified in the configuration file or you wish to overwrite any of the references.
    --fasta                       Path to Fasta reference
    --gff                         Path to GFF reference
    --genus                       Genus information for use in CheckM

    External databases:
    --nt_db                       NCBI Nt database (BLAST)
    --uniprot_db                  Uniprot proteomes database (diamond) !!! time consuming !!!
    --uniprot_taxids              Sequence id to taxa id mapping file
    --kraken_db                   Kraken database
    --eggnog_db                   EggNOG v4.5.1 database for emapper-1.0.3
    --kofam_profile               KOfam profile database
    --kofam_kolist                KOfam ko_list file
    --augustus_species            Augustus species, default 'saccharomyces'
    --eukcc_db                    EukCC database

    Trimming options:
    --notrim                      Specifying --notrim will skip the adapter trimming step.
    --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.
    --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
    --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
    --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
    --three_prime_clip_r2 [int]   Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed

    Mapping options:
    --allow_multi_align           Secondary alignments and unmapped reads are also reported in addition to primary alignments
    --saveAlignedIntermediates    Save the intermediate BAM files from the Alignment step  - not done by default

    Assembly options:
    --no_normalize                Specifying --no_normalize will skip the reads normalizing step.

    Quast options:
    --euk                         Euk genome
    --fungus                      Fungal genome

    Taxa annotation options:
    --evalue                      E-value for blasting NCBI-nt and uniprot reference proteomes database (default=1e-25)

    Diamond options:
    --blockSize                   Sequence block size in billions of letters (default=2.0)

    ARG related options:
    --acquired                    Enable ARG analysis
    --point                       Enable point mutation analysis
    --only_known                  Only analyze known SNPs
    --resfinder_db                Database path for resfinder
    --pointfinder_db              Database path for pointfinder
    --pointfinder_species         Species for pointfinder, default 'escherichia_coli'

    Output options:
    --outdir                      The output directory where the results will be saved
    --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
    --maxMultiqcEmailFileSize     Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)

    AWSBatch options:
    --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
    --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// default values
params.genome = false
params.nanopore = false
params.single_end = false
params.fasta = false
params.gff = false
params.vcf = false
params.notrim = false
params.saveTrimmed = false
params.allow_multi_align = false
params.saveAlignedIntermediates = false
params.no_normalize = false
params.euk = false
params.fungus = false
params.remap = false
params.genus = null
params.nt_db = null
params.kraken_db = null
params.readPaths = null
params.uniprot_db = null
params.uniprot_taxids = null
params.eggnog_db = null
params.eukcc_db = null
params.snv = false
params.cnv = false
params.saturation = false
params.bulk = false
params.ass = false
params.evalue = 1e-25
params.blockSize = 2.0
params.acquired = false
params.point = false
params.only_known = true
params.pointfinder_species = "escherichia_coli"
params.resfinder_db = null
params.pointfinder_db = null
params.kofam_profile = null
params.kofam_kolist = null
params.augustus_species = "saccharomyces"
params.split = false
params.split_bac_level = null
params.split_euk_level = null

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Configurable reference genomes
fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if ( params.fasta ){
    fasta = file(params.fasta)
    ref = params.fasta - ~/(\.fasta)?(\.fna)?(\.fa)?$/
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//

gff = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
if ( params.gff ) {
    gff = file(params.gff)
    if( !gff.exists() ) exit 1, "GFF file not found: ${params.gff}"
}

graph_vcf = false
if ( params.vcf ) {
    graph_vcf = file(params.vcf)
    if ( !graph_vcf.exists()) exit 1, "VCF file to construct graph not found: ${params.graph_vcf}"
}

single_end = false
if ( params.nanopore ) {
    single_end = true
} else {
    single_end = params.single_end
}

euk = false
if ( params.fungus || params.euk ) {
    euk = true
}

// Configurable nt database
nt_db = false
if ( params.nt_db ) {
    nt_db = file(params.nt_db)
    if( !nt_db.exists() ) exit 1, "NT database not found: ${params.nt_db}"
}

// Configurable uniprot proteomes database
uniprot_db = false
if ( params.uniprot_db ) {
    uniprot_db = file(params.uniprot_db)
    if ( !uniprot_db.exists() ) exit 1, "Uniprot proteomes database not found: ${params.uniprot_db}"
} else {
    uniprot_db = file("/dev/null")
}

//uniprot_taxids
uniprot_taxids = false
if ( params.uniprot_taxids ) {
    uniprot_taxids = file(params.uniprot_taxids)
    if ( !uniprot_taxids.exists() ) exit 1, "Uniprot proteomes seq2tax mapping file not found: ${params.uniprot_taxids}"
} else {
    uniprot_taxids = file("/dev/null")
}

// Configurable kraken database
kraken_db = false
if ( params.kraken_db ) {
    kraken_db = file(params.kraken_db)
    if( !kraken_db.exists() ) exit 1, "Kraken database not found: ${params.kraken_db}"
}

// Configurable eggNOG database
eggnog_db = false
if ( params.eggnog_db ) {
    eggnog_db = file(params.eggnog_db)
    if( !eggnog_db.exists() ) exit 1, "EggNOG database not found: ${params.eggnog_db}"
}

// Configure resfinder_db
resfinder_db = false
if ( params.resfinder_db ) {
    resfinder_db = file(params.resfinder_db)
    if( !resfinder_db.exists() ) exit 1, "ResFinder database not found: ${params.resfinder_db}"
}

// Configure pointfinder_db
pointfinder_db = false
if ( params.pointfinder_db ) {
    pointfinder_db = file(params.pointfinder_db)
    if( !pointfinder_db.exists() ) exit 1, "PointFinder database not found: ${params.pointfinder_db}"
}

// Configure EukCC database
eukcc_db = false
if ( params.eukcc_db ) {
    eukcc_db  = file(params.eukcc_db)
    if ( !eukcc_db.exists() ) exit 1, "EukCC database not found: ${params.eukcc_db}"
}

// Configure KOfam search database
kofam_profile = false
if ( params.kofam_profile ) {
    kofam_profile = file(params.kofam_profile)
    if( !kofam_profile.exists() ) exit 1, "KOfam profile database not found: ${params.kofam_profile}"
}

kofam_kolist = false
if ( params.kofam_kolist ) {
    kofam_kolist = file(params.kofam_kolist)
    if( !kofam_kolist.exists() ) exit 1, "KOfam ko_list file not found: ${params.kofam_kolist}"
}

custom_runName = workflow.runName

if( workflow.profile == 'awsbatch') {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
    // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Stage config files
ch_multiqc_config1 = Channel.fromPath(params.multiqc_config)
ch_multiqc_config2 = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

// Custom trimming options
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// mode
denovo = (params.genome && params.genomes[ params.genome ].bowtie2) || params.fasta ? false : true
bowtie2 = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false


/*
 * Create a channel for input read files
 */
if(params.readPaths){
    if(single_end){
        read_files_fastqc = read_files_trimming =
        Channel.from(params.readPaths, checkIfExists: true)
            .map { row -> def meta=[:];
                    meta.id = row[0].replaceFirst(~/\.[^\.]+$/, '');
                    meta.single_end = single_end;
                    [meta, [file(row[1][0]), file(row[1][1])]]}
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
    } else {
        read_files_fastqc = read_files_trimming =
        Channel.from(params.readPaths)
            .map { row -> def meta=[:];
                    meta.id = row[0].replaceFirst(~/\.[^\.]+$/, '');
                    meta.single_end = single_end;
                    [meta, [file(row[1][0]), file(row[1][1])]]}
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
    }
} else {
    if (single_end) {
        read_files_fastqc = read_files_trimming =
        Channel.fromFilePairs( params.reads, size:1, checkIfExists: true)
            .map { it ->
                def meta = [:];
                meta.id = it[0].replaceFirst(~/\.[^\.]+$/, '');
                meta.single_end = single_end;
                [meta, [file(it[1][0])]]}

        } else {
            read_files_fastqc = read_files_trimming =
            Channel.fromFilePairs( params.reads, size:2, checkIfExists: true)
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
summary['Fasta Ref']        = params.fasta
summary['Data Type']        = single_end ? 'Single-End' : 'Paired-End'
summary['Bulk']             = params.bulk ? 'Yes' : 'No'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if( params.notrim ){
    summary['Trimming Step'] = 'Skipped'
} else {
    summary["Trim 5' R1"] = params.clip_r1
    summary["Trim 5' R2"] = params.clip_r2
    summary["Trim 3' R1"] = params.three_prime_clip_r1
    summary["Trim 3' R2"] = params.three_prime_clip_r2
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

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-scgs-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'gongyh/nf-core-scgs Workflow Summary'
    section_href: 'https://github.com/gongyh/nf-core-scgs'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v != null ? v : '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

    return yaml_file
}

// Import modules from nf-core
include { FASTQC                } from './modules/nf-core/fastqc/main'
include { TRIMGALORE            } from './modules/nf-core/trimgalore/main'
include { BOWTIE2_BUILD         } from './modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN         } from './modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN        } from './modules/nf-core/minimap2/align/main'
include { VG_CONSTRUCT          } from './modules/nf-core/vg/construct/main'
include { VG_INDEX              } from './modules/nf-core/vg/index/main'
include { QUALIMAP_BAMQC        } from './modules/nf-core/qualimap/bamqc/main'

// Import modules from local
// include { GET_SOFTWARE_VERSIONS } from './modules/local/get_software_versions'
include { SAVE_REFERENCE        } from './modules/local/save_reference'
include { KRAKEN                } from './modules/local/kraken'
include { SATURATION            } from './modules/local/saturation'
include { VG_CALL               } from './modules/local/vg_call'
include { SAMTOOLS              } from './modules/local/samtools'
include { PRESEQ                } from './modules/local/preseq'
include { INDELREALIGN          } from './modules/local/indelrealign'
include { MONOVAR               } from './modules/local/monovar'
include { ANEUFINDER            } from './modules/local/aneufinder'
include { CIRCLIZE              } from './modules/local/circlize'
include { NORMALIZE             } from './modules/local/normalize'
include { CANU                  } from './modules/local/canu'
include { SPADES                } from './modules/local/spades'
include { BOWTIE2_REMAP         } from './modules/local/bowtie2_remap'
include { REMAP                 } from './modules/local/remap'
include { QUAST_REF             } from './modules/local/quast_ref'
include { QUAST_DENOVO          } from './modules/local/quast_denovo'
include { CHECKM_LINEAGEWF      } from './modules/local/checkm/lineagewf/main'
include { BLASTN                } from './modules/local/blastn'
include { DIAMOND_BLASTX        } from './modules/local/diamond/blastx/main'
include { BLOBTOOLS             } from './modules/local/blobtools'
include { REBLOBTOOLS           } from './modules/local/reblobtools'
include { ACDC                  } from './modules/local/acdc'
include { TSNE                  } from './modules/local/tsne'
include { PROKKA                } from './modules/local/prokka/main'
include { PRODIGAL              } from './modules/local/prodigal/main'
include { AUGUSTUS              } from './modules/local/augustus'
include { EUKCC                 } from './modules/local/eukcc'
include { MULTIQC_DENOVO        } from './modules/local/multiqc/multiqc_denovo'
include { MULTIQC_REF           } from './modules/local/multiqc/multiqc_ref'
include { EGGNOG                } from './modules/local/eggnog'
include { KOFAM                 } from './modules/local/kofam'
include { RESFINDER             } from './modules/local/resfinder'
include { POINTFINDER           } from './modules/local/pointfinder'
include { SPLIT_CHECKM          } from './modules/local/split/split_checkm'
include { SPLIT_CHECKM_EUKCC    } from './modules/local/split/split_checkm_eukcc'
include { OUTPUT_DOCUMENTATION  } from './modules/local/output_documentation'


workflow {
    ch_versions = Channel.empty()

    // FASTQC
    FASTQC(read_files_fastqc)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first().ifEmpty(null))

    if (params.fasta && params.gff) {
        SAVE_REFERENCE(fasta, gff)
    }

    if (bowtie2) {
        bowtie2_index = [bowtie2, file(bowtie2)]
    } else {
        if (params.fasta) {
            BOWTIE2_BUILD([ref, fasta])
            bowtie2_index = BOWTIE2_BUILD.out.index
        }
    }
    // TRIM_GALORE
    trimmed_reads = Channel.empty()
    if (params.notrim) {
        trimmed_reads = read_files_trimming.map {name, reads -> reads}
        trimgalore_results1 = file('/dev/null')
        trimgalore_results2 = file('/dev/null')
        trimgalore_fastqc_reports1 = file('/dev/null')
        trimgalore_fastqc_reports2 = file('/dev/null')
    } else {
        TRIMGALORE(read_files_trimming)
        trimmed_reads = TRIMGALORE.out.reads
    }
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first().ifEmpty(null))

    // KRAKEN
    if (params.kraken_db) {
        KRAKEN(trimmed_reads, kraken_db)
    }
    SATURATION(trimmed_reads)

    // ALIGN
    if (denovo == false) {
        bb_bam = Channel.empty()
        if (params.nanopore) {
            MINIMAP2_ALIGN(trimmed_reads,
                            fasta,
                            true,
                            false,
                            false)
            MINIMAP2_ALIGN.out.bam
                        .set { bb_bam }
        } else {
            BOWTIE2_ALIGN(trimmed_reads,
                            bowtie2_index,
                            false,
                            true)
            BOWTIE2_ALIGN.out.bam
                        .set { bb_bam }
        }
    }

    // VG


    if (params.fasta && params.gff) {
        SAMTOOLS(bb_bam, SAVE_REFERENCE.out.bed)
        PRESEQ(SAMTOOLS.out.bed)
        QUALIMAP_BAMQC(SAMTOOLS.out.bam , gff)
        INDELREALIGN(SAMTOOLS.out.bam, fasta)
        MONOVAR(INDELREALIGN.out.bam.map { meta,bam -> return bam}.collect(), INDELREALIGN.out.bai.map {meta, bai -> return bai}.collect(), fasta)
        ANEUFINDER(SAMTOOLS.out.bam.map { meta,bam -> return bam}.collect(), SAMTOOLS.out.bai.map {meta, bai -> return bai}.collect())
        CIRCLIZE(SAMTOOLS.out.bed, SAVE_REFERENCE.out.bed)
    }

    // NORMALIZE
    ctg200 = Channel.empty()
    ctg = Channel.empty()
    if (params.no_normalize) {
        trimmed_reads.set { normalized_reads }
    } else {
        if (params.ass) {
            NORMALIZE(trimmed_reads)
            if (params.nanopore) {
                CANU(NORMALIZE.out.reads)
                ctg200 = CANU.out.ctg200
                ctg = CANU.out.ctg
            } else {
                SPADES(NORMALIZE.out.reads)
                ctg200 = SPADES.out.ctg200
                ctg = SPADES.out.ctg
            }
        }
    }

    // REMAP
    if (params.remap) {
        BOWTIE2_REMAP(ctg200)
        REMAP(trimmed_reads, BOWTIE2_REMAP.out.index.collet())
    }

    // QUAST
    quast_denovo = Channel.empty()
    if (denovo == false) {
        QUAST_REF(fasta, gff, ctg, SAMTOOLS.out.bam.collect(), SAMTOOLS.out.bai.collect(), euk)
    } else {
        QUAST_DENOVO(ctg, euk)
        quast_denovo = QUAST_DENOVO.out.report
    }

    /**
    // CHECKM
    if (!params.euk) {
        CHECKM_LINEAGEWF(ctg.collect())
    }

    // NT
    if (params.nt_db) {
        BLASTN(ctg200, nt_db)
        DIAMOND_BLASTX(BLASTN.out.contigs, BLASTN.out.nt, uniprot_db, uniprot_taxids)
        BLOBTOOLS(DIAMOND_BLASTX.out.contigs, DIAMOND_BLASTX.out.nt, DIAMOND_BLASTX.out.real, DIAMOND_BLASTX.out.uniprot)
        if (params.remap) {
            REBOLBTOOLS(DIAMOND_BLASTX.out.contigs, DIAMOND_BLASTX.out.nt, DIAMOND_BLASTX.out.real, DIAMOND_BLASTX.out.uniprot, REMAP.out.bam.collect())
        }
        ACDC(BLOBTOOLS.out.contigs, kraken_db, BLOBTOOLS.out.tax)
    }
    TSNE(ctg)

    faa = Channel.empty()
    if (!params.euk) {
        PROKKA(ctg)
        PRODIGAL(ctg)
        faa = PROKKA.out.faa
        prokka_for_mqc1 = PROKKA.out.prokka_for_split
        prokka_for_mqc2 = PROKKA.out.prokka_for_split
        prokka_for_split = PROKKA.out.prokka_for_split
    } else {
        prokka_for_mqc1 = file('/dev/null')
        prokka_for_mqc2 = file('/dev/null')
        prokka_for_split = file('/dev/null')
        AUGUSTUS(ctg)
        faa = AUGUSTUS.out.faa
        EUKCC(faa, eukcc_db)
    }

    ch_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_versions }

    GET_SOFTWARE_VERSIONS (
        ch_versions.map { it }.collect()
    )

    // MULTIQC
    preseq_for_multiqc = file('/dev/null')

    MULTIQC_REF(
                ch_multiqc_config1,
                FASTQC.out.zip.collect{it[1]}.ifEmpty([]),
                TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]),
                GET_SOFTWARE_VERSIONS.out.yaml.collect(),
                SAMTOOLS.out.stats.collect(),
                preseq_for_multiqc.collect(),
                QUALIMAP_BAMQC.out.results,
                QUAST_REF.out.report,
                prokka_for_mqc1.collect(),
                KRAKEN.out.report.collect(),
                create_workflow_summary(summary)
                )
    MULTIQC_DENOVO(
                ch_multiqc_config2,
                FASTQC.out.zip.collect{it[1]}.ifEmpty([]),
                TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]),
                GET_SOFTWARE_VERSIONS.out.yaml.collect(),
                quast_denovo.ifEmpty('/dev/null'),
                prokka_for_mqc2.collect(),
                KRAKEN.out.report.collect(),
                create_workflow_summary(summary)
                )

    if (params.eggnog) {
        EGGNOG(faa, eggnog_db)
    }

    if (params.kofam_profile && params.kofam_kolist) {
        KOFAM(faa, kofam_profile, kofam_kolist)
    }

    if (!params.euk && params.acquired && resfinder_db) {
        RESFINDER(ctg, resfinder_db)
    }

    if (!params.euk && params.point && pointfinder_db) {
        POINTFINDER(ctg, pointfinder_db)
    }

    if (params.eukcc_db) {
        SPLIT_CHECKM_EUKCC(
            ctg200.collect(),
            BLOBTOOLS.out.tax.collect(),
            PROKKA.out.prokka_for_split.collect(),
            KOFAM.out.txt.collect(),
            eukcc_db,
            params.split_bac_level ? params.split_bac_level : "genus",
            params.split_euk_level ? params.split_euk_level : "genus"
        )
    } else {
        SPLIT_CHECKM(
            ctg200.collect(),
            BLOBTOOLS.out.tax.collect(),
            PROKKA.out.prokka_for_split.collect(),
            KOFAM.out.txt.collect()
        )
    }

    OUTPUT_DOCUMENTATION(ch_output_docs)
    */

}

/*
 * Completion e-mail notification
 */
workflow.onComplete {
    //NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    //NfcoreTemplate.summary(workflow, params, log)
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
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  gongyh/nf-core-scgs v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}