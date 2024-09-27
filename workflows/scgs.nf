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
    --bbmap                       Enable bbmap to remove host-derived contamination
    --doubletd                    Enable detection of doublet
    --remap                       Remap trimmed reads to contigs
    --acdc                        Enable acdc
    --saturation                  Enable sequencing saturation analysis
    --ass                         Assemble using SPAdes
    --genomad                     Enable genomad analysis
    --blastn                      Enable NCBI Nt database annotation
    --blob                        Enable Blobtools analysis
    --kraken                      Enable Kraken2 annotation
    --eggnog                      Enable EggNOG database annotation
    --kofam                       Enable KEGG Ortholog annotation
    --checkm2                     Enable CheckM2 analysis
    --gtdbtk                      Enable gtdbtk analysis
    --split                       Split the draft genomes and annotation(Bacteria)
    --split_euk                   Split the draft genomes and annotation(Eukaryota)
    --split_bac_level             Level of split for Bacteria
    --split_euk_level             Level of split for Eukaryota
    --graphbin                    Enable graphbin to bin
    --pangenome                   Enable pangenome analysis
    --completeness                Calculate the completeness of assembling contigs using pan-genomic methods based on core genes
    --tree                        Draw a phylogenetic tree

    References:                   If not specified in the configuration file or you wish to overwrite any of the references.
    --fasta                       Path to Fasta reference
    --gff                         Path to GFF reference
    --genus                       Genus information for use in CheckM

    External databases:
    --genomad_db                  geNomad database
    --prokka_proteins             FASTA file of trusted proteins to first annotate from (optional)
    --nt_db                       NCBI Nt database (BLAST)
    --blob_db                     Blobtools nodesDB.txt
    --uniprot_db                  Uniprot proteomes database (diamond) !!! time consuming !!!
    --uniprot_taxids              Sequence id to taxa id mapping file
    --kraken_db                   Kraken2 database
    --eggnog_db                   EggNOG v4.5.1 database for emapper-1.0.3
    --kofam_profile               KOfam profile database
    --kofam_kolist                KOfam ko_list file
    --augustus_species            Augustus species, default 'saccharomyces'
    --eukcc_db                    EukCC database
    --checkm2_db                  CheckM2 database
    --mgpg_db                     Microbiome graph pangenome database
    --gtdb                        GTDB database
    --ref                         Specify the reference sequence for bbmap

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

    Pangenome options:
    --genusName                   Genus Name
    --coreGenesFile               Core genes txt file

    Taxa annotation options:
    --evalue                      E-value for blasting NCBI-nt and uniprot reference proteomes database (default=1e-25)

    Diamond options:
    --blockSize                   Sequence block size in billions of letters (default=2.0)

    ARG related options:
    --acquired                    Enable ARG analysis
    --point                       Enable point mutation analysis
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
params.acdc = false
params.snv = false
params.cnv = false
params.bbmap = false
params.doubletd = false
params.saturation = false
params.bulk = false
params.genomad = false
params.ass = false
params.kofam = true
params.blastn = true
params.blob = true
params.kraken = true
params.eggnog = true
params.checkm2 = true
params.gtdbtk = true
params.acquired = false
params.point = false
params.split = false
params.split_euk = false
params.graphbin = false
params.pangenome = false
params.completeness = false
params.tree = false
params.genusName = null
params.coreGenesFile = null
params.genus = null
params.genomad_db = null
params.prokka_proteins = null
params.nt_db = null
params.blob_db = null
params.kraken_db = null
params.kofam_profile = null
params.kofam_kolist = null
params.readPaths = null
params.uniprot_db = null
params.uniprot_taxids = null
params.eggnog_db = null
params.eukcc_db = null
params.checkm2_db = null
params.gtdb = null
params.ref = null
params.evalue = 1e-25
params.blockSize = 2.0
params.split_bac_level = "genus"
params.split_euk_level = "genus"
params.pointfinder_species = "escherichia_coli"
params.augustus_species = "saccharomyces"


// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Configurable reference genomes
fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    fasta = file(params.fasta)
    ref = params.fasta - ~/(\.fasta)?(\.fna)?(\.fa)?$/
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}

gff = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
if (params.gff) {
    gff = file(params.gff)
    if( !gff.exists() ) exit 1, "GFF file not found: ${params.gff}"
}

graph_vcf = false
if (params.vcf) {
    graph_vcf = file(params.vcf)
    if ( !graph_vcf.exists()) exit 1, "VCF file to construct graph not found: ${params.graph_vcf}"
} else {
    graph_vcf = file("/dev/null")
}

single_end = false
if (params.nanopore) {
    single_end = true
} else {
    single_end = params.single_end
}

euk = false
if (params.fungus || params.euk) {
    euk = true
}

// Configurable genomad database
genomad_db = false
if (params.genomad_db) {
    genomad_db = file(params.genomad_db)
    if( !genomad_db.exists() ) exit 1, "Genomad database not found: ${params.genomad_db}"
} else {
    genomad_db = file("/dev/null")
}

// Prokka trusted proteins database
prokka_proteins = []
if (params.prokka_proteins) {
    faa = file(params.prokka_proteins)
    if( !prokka_proteins.exists() ) exit 1, "Protein database not found: ${params.prokka_proteins}"
    prokka_proteins = [faa]
}

// Configurable nt database
nt_db = false
if (params.nt_db) {
    nt_db = file(params.nt_db)
    if( !nt_db.exists() ) exit 1, "NT database not found: ${params.nt_db}"
} else {
    nt_db = file("/dev/null")
}

// Configurable uniprot proteomes database
uniprot_db = false
if (params.uniprot_db) {
    uniprot_db = file(params.uniprot_db)
    if ( !uniprot_db.exists() ) exit 1, "Uniprot proteomes database not found: ${params.uniprot_db}"
} else {
    uniprot_db = file("/dev/null")
}

//uniprot_taxids
uniprot_taxids = false
if (params.uniprot_taxids) {
    uniprot_taxids = file(params.uniprot_taxids)
    if ( !uniprot_taxids.exists() ) exit 1, "Uniprot proteomes seq2tax mapping file not found: ${params.uniprot_taxids}"
} else {
    uniprot_taxids = file("/dev/null")
}

// Configurable kraken database
kraken_db = false
if (params.kraken_db) {
    kraken_db = file(params.kraken_db)
    if( !kraken_db.exists() ) exit 1, "Kraken database not found: ${params.kraken_db}"
} else {
    kraken_db = file("/dev/null")
}

// Configurable Blobtools nodesDB.txt
blob_db = false
if (params.blob_db) {
    blob_db = file(params.blob_db)
    if( !blob_db.exists() ) exit 1, "Blobtools nodesDB.txt not found: ${params.blob_db}"
} else {
    blob_db = file("/dev/null")
}

// Configurable eggNOG database
eggnog_db = false
if (params.eggnog_db) {
    eggnog_db = file(params.eggnog_db)
    if( !eggnog_db.exists() ) exit 1, "EggNOG database not found: ${params.eggnog_db}"
} else {
    eggnog_db = file("/dev/null")
}

// Configure EukCC database
eukcc_db = false
if (params.eukcc_db) {
    eukcc_db  = file(params.eukcc_db)
    if ( !eukcc_db.exists() ) exit 1, "EukCC database not found: ${params.eukcc_db}"
} else {
    eukcc_db = file("/dev/null")
}

// Configure Checkm2 database
checkm2_db = false
if (params.checkm2_db) {
    checkm2_db  = file(params.checkm2_db)
    if ( !checkm2_db.exists() ) exit 1, "CheckM2 database not found: ${params.checkm2_db}"
} else {
    checkm2_db = file("/dev/null")
}

// Configure GTDB database
gtdb = false
if (params.gtdb) {
    gtdb  = file(params.gtdb)
    if ( !gtdb.exists() ) exit 1, "GTDB database not found: ${params.gtdb}"
} else {
    gtdb = file("/dev/null")
}

// Configure pangenome database
mgpg_db = false
if (params.mgpg_db) {
    mgpg_db  = file(params.mgpg_db)
    if ( !mgpg_db.exists() ) exit 1, "Graph pangenome database not found: ${params.mgpg_db}"
} else {
    mgpg_db = file("/dev/null")
}

// Configure core genes file
coreGenesFile = false
if (params.coreGenesFile) {
    coreGenesFile  = file(params.coreGenesFile)
    if ( !coreGenesFile.exists() ) exit 1, "Core genes file not found: ${params.coreGenesFile}"
} else {
    coreGenesFile = file("/dev/null")
}

// Configure reference sequence
ref = false
if (params.ref) {
    ref  = file(params.ref)
    if ( !ref.exists() ) exit 1, "Ref not found: ${params.ref}"
} else {
    ref = file("/dev/null")
}

// Configure KOfam search database
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

custom_runName = workflow.runName

if(workflow.profile == 'awsbatch') {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
    // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
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
include { FASTQC                } from '../modules/nf-core/fastqc/main'
include { BOWTIE2_BUILD         } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN         } from '../modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN        } from '../modules/nf-core/minimap2/align/main'
include { QUALIMAP_BAMQC        } from '../modules/nf-core/qualimap/bamqc/main'
include { GENOMAD_ENDTOEND      } from '../modules/nf-core/genomad/endtoend/main'
include { MULTIQC               } from '../modules/nf-core/multiqc/main'


// Import modules from local
include { SAVE_REFERENCE        } from '../modules/local/save_reference'
include { TRIMGALORE            } from '../modules/local/trimgalore'
include { KTUPDATETAXONOMY      } from '../modules/local/ktupdatetaxonomy'
include { KRAKEN                } from '../modules/local/kraken'
include { SATURATION            } from '../modules/local/saturation'
include { SAMTOOLS              } from '../modules/local/samtools'
include { PRESEQ                } from '../modules/local/preseq'
include { GTDBTK                } from '../modules/local/gtdbtk'
include { BBMAP_ALIGN           } from '../modules/local/bbmap_align'
include { INDELREALIGN          } from '../modules/local/indelrealign'
include { MONOVAR               } from '../modules/local/monovar'
include { DOUBLETD              } from '../modules/local/doubletd'
include { ANEUFINDER            } from '../modules/local/aneufinder'
include { CIRCLIZE              } from '../modules/local/circlize'
include { NORMALIZE             } from '../modules/local/normalize'
include { BBNORM                } from '../modules/local/bbnorm'
include { CANU                  } from '../modules/local/canu'
include { SPADES                } from '../modules/local/spades'
include { COMPLETENESS          } from '../modules/local/pangenome/completeness'
include { TREE                  } from '../modules/local/pangenome/tree'
include { QUAST_REF             } from '../modules/local/quast_ref'
include { QUAST_DENOVO          } from '../modules/local/quast_denovo'
include { BOWTIE2_REMAP         } from '../modules/local/bowtie2_remap'
include { REMAP                 } from '../modules/local/remap'
include { CHECKM_LINEAGEWF      } from '../modules/local/checkm_lineagewf'
include { CHECKM2               } from '../modules/local/checkm2'
include { BLASTN                } from '../modules/local/blastn'
include { DIAMOND_BLASTX        } from '../modules/local/diamond_blastx'
include { BLOBTOOLS             } from '../modules/local/blobtools'
include { REBLOBTOOLS           } from '../modules/local/reblobtools'
include { ACDC                  } from '../modules/local/acdc'
include { TSNE                  } from '../modules/local/tsne'
include { PROKKA                } from '../modules/local/prokka'
include { PRODIGAL              } from '../modules/local/prodigal'
include { METARON               } from '../modules/local/metaron'
include { AUGUSTUS              } from '../modules/local/augustus'
include { EUKCC                 } from '../modules/local/eukcc'
include { EGGNOG                } from '../modules/local/eggnog'
include { KOFAMSCAN             } from '../modules/local/kofamscan'
include { STARAMR               } from '../modules/local/staramr'
include { SPLIT_CHECKM          } from '../modules/local/split_checkm'
include { SPLIT_CHECKM_EUKCC    } from '../modules/local/split_checkm_eukcc'
include { GRAPHBIN              } from '../modules/local/graphbin'
include { OUTPUT_DOCUMENTATION  } from '../modules/local/output_documentation'
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions/main'


/** subworkflow */
include { VG                    } from '../subworkflows/local/vg'

// MULTIQC
def multiqc_report = []

workflow SCGS {
    ch_versions = Channel.empty()

    // FASTQC
    ch_multiqc_fastqc = Channel.empty()
    FASTQC ( read_files_fastqc )
    ch_versions       = ch_versions.mix(FASTQC.out.versions)
    ch_multiqc_fastqc = FASTQC.out.zip

    // SAVE_REFERENCE
    if ( params.fasta ) {
        if ( params.gff ) {
            SAVE_REFERENCE ( fasta, gff )
        } else {
            SAVE_REFERENCE ( fasta, file("/dev/null") )
        }
    }

    if (bowtie2) {
        bowtie2_index = [bowtie2, file(bowtie2)]
    } else {
        if (params.fasta) {
            BOWTIE2_BUILD ( [ref, fasta] )
            bowtie2_index = BOWTIE2_BUILD.out.index
        }
    }
    // TRIM_GALORE
    trimmed_reads = Channel.empty()
    ch_multiqc_trim_log = Channel.empty()
    ch_multiqc_trim_zip = Channel.empty()
    if (params.notrim) {
        if (params.bbmap) {
            BBMAP_ALIGN (
                read_files_trimming.map{name, reads -> reads},
                ref
            )
            ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)
            trimmed_reads = BBMAP_ALIGN.out.clean_fastq
        } else {
            trimmed_reads = read_files_trimming.map{name, reads -> reads}
        }
    } else {
        TRIMGALORE ( read_files_trimming )
        ch_multiqc_trim_log = TRIMGALORE.out.log
        ch_multiqc_trim_zip = TRIMGALORE.out.zip
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
        if (params.bbmap) {
            BBMAP_ALIGN (
                TRIMGALORE.out.reads,
                ref
            )
            ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)
            trimmed_reads = BBMAP_ALIGN.out.clean_fastq
        } else {
            trimmed_reads = TRIMGALORE.out.reads
        }
    }

    // KRAKEN
    ch_multiqc_kraken = Channel.empty()
    if (params.kraken) {
        KTUPDATETAXONOMY ()
        KRAKEN (
            trimmed_reads,
            kraken_db,
            KTUPDATETAXONOMY.out.taxonomy
        )
        ch_multiqc_kraken = KRAKEN.out.report
    }

    // SATURATION
    if (params.saturation) {
        SATURATION ( trimmed_reads )
    }

    // ALIGN
    if (denovo == false) {
        bb_bam = Channel.empty()
        if ( params.nanopore ) {
            MINIMAP2_ALIGN (
                trimmed_reads,
                fasta,
                true,
                false,
                false
            )
            MINIMAP2_ALIGN.out.bam.set{ bb_bam }
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
        } else {
            BOWTIE2_ALIGN (
                trimmed_reads,
                bowtie2_index,
                false,
                true
            )
            BOWTIE2_ALIGN.out.bam.set{ bb_bam }
            ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
        }
    }

    // VG
    if ( params.fasta && params.vcf ) {
        VG (
            fasta,
            trimmed_reads,
            graph_vcf
        )
        ch_versions = ch_versions.mix(VG.out.ch_versions)
    }

    ch_multiqc_samtools = Channel.empty()
    ch_multiqc_preseq   = Channel.empty()
    ch_multiqc_qualimap = Channel.empty()
    quast_bam = Channel.empty()
    quast_bai = Channel.empty()
    if ( params.fasta ) {
        SAMTOOLS (
            bb_bam,
            SAVE_REFERENCE.out.bed
        )
        SAMTOOLS.out.bam.set{quast_bam}
        SAMTOOLS.out.bai.set{quast_bai}
        ch_versions = ch_versions.mix(SAMTOOLS.out.versions)
        ch_multiqc_samtools = SAMTOOLS.out.stats
        if (!params.nanopore) {
            PRESEQ ( SAMTOOLS.out.bed )
            ch_versions = ch_versions.mix(PRESEQ.out.versions)
            ch_multiqc_preseq = PRESEQ.out.txt
        }
        if ( params.gff ) {
            QUALIMAP_BAMQC (
                SAMTOOLS.out.bam,
                gff
            )
            ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)
            ch_multiqc_qualimap = QUALIMAP_BAMQC.out.results
        }
        if (params.snv && !params.nanopore) {
            INDELREALIGN (
                SAMTOOLS.out.bam,
                fasta
            )
            ch_versions = ch_versions.mix(INDELREALIGN.out.versions)
        }
        if (!params.bulk && params.snv && !params.nanopore) {
            MONOVAR (
                INDELREALIGN.out.bam.collect{it[1]},
                INDELREALIGN.out.bai.collect{it[1]},
                fasta
            )
            ch_versions = ch_versions.mix(MONOVAR.out.versions)
            if ( params.doubletd ) {
                DOUBLETD ( MONOVAR.out.vcf )
            }
        }
        if (!params.bulk && params.cnv && !single_end && !params.nanopore) {
            ANEUFINDER (
                SAMTOOLS.out.bam.collect{it[1]},
                SAMTOOLS.out.bai.collect{it[1]}
            )
            ch_versions = ch_versions.mix(ANEUFINDER.out.versions)
        }
        CIRCLIZE (
            SAMTOOLS.out.bed,
            SAVE_REFERENCE.out.bed
        )
        ch_versions = ch_versions.mix(CIRCLIZE.out.versions)
    }

    // ASSEMBLY
    ctg200 = Channel.empty()
    ctg = Channel.empty()
    if ( params.ass ) {
        // NORMALIZE
        if ( params.no_normalize ) {
            trimmed_reads.set{ normalized_reads }
        } else {
            /**
            NORMALIZE(trimmed_reads)
            normalized_reads = NORMALIZE.out.reads
            */
            BBNORM(trimmed_reads)
            normalized_reads = BBNORM.out.fastq
            ch_versions = ch_versions.mix(BBNORM.out.versions)
        }
        contig = Channel.empty()
        contig_path = Channel.empty()
        contig_graph = Channel.empty()
        if (params.nanopore) {
            CANU(normalized_reads)
            ctg200 = CANU.out.ctg200
            ctg = CANU.out.ctg
            ch_versions = ch_versions.mix(CANU.out.versions)
        } else {
            SPADES(normalized_reads)
            contig = SPADES.out.contig
            contig_path = SPADES.out.contig_path
            contig_graph = SPADES.out.contig_graph
            ctg200 = SPADES.out.ctg200
            ctg = SPADES.out.ctg
            ch_versions = ch_versions.mix(SPADES.out.versions)
        }
    }

    // GENOMAD
    if ( params.genomad ) {
        GENOMAD_ENDTOEND(
          ctg,
          genomad_db
        )
        ch_versions = ch_versions.mix(GENOMAD_ENDTOEND.out.versions)
    }

    // REMAP
    if (params.remap && !params.nanopore) {
        BOWTIE2_REMAP(ctg200)
        REMAP (
            trimmed_reads,
            BOWTIE2_REMAP.out.index.collect{it[1]},
            params.allow_multi_align
        )
    }

    // QUAST
    ch_multiqc_quast = Channel.empty()
    if (denovo == false) {
        QUAST_REF (
            fasta,
            gff,
            ctg.collect{it[1]},
            quast_bam.collect{it[1]},
            quast_bai.collect{it[1]},
            euk,
            params.fungus
        )
        ch_multiqc_quast = QUAST_REF.out.tsv
        ch_versions = ch_versions.mix(QUAST_REF.out.versions)
    } else {
        QUAST_DENOVO (
            ctg.collect{it[1]},
            euk,
            params.fungus
        )
        ch_multiqc_quast = QUAST_DENOVO.out.tsv
        ch_versions = ch_versions.mix(QUAST_DENOVO.out.versions)
    }

    // CHECKM_LINEAGEWF
    ch_multiqc_checkm = Channel.empty()
    if (!euk) {
        CHECKM_LINEAGEWF (
            ctg.collect{it[1]},
            params.genus ? true : false
        )
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions)
        ch_multiqc_checkm = CHECKM_LINEAGEWF.out.mqc_tsv
    }

    // CHECKM2
    ch_multiqc_checkm2 = Channel.empty()
    if (!euk && params.checkm2) {
        CHECKM2 (
            ctg.collect{it[1]},
            checkm2_db
        )
        ch_versions = ch_versions.mix(CHECKM2.out.versions)
        ch_multiqc_checkm2 = CHECKM2.out.mqc_tsv
    }

    if (params.blastn) {
        // BLASTN
        BLASTN (
            ctg200,
            nt_db,
            params.evalue
        )
        ch_versions = ch_versions.mix(BLASTN.out.versions)

        // DIAMOND_BLASTS
        DIAMOND_BLASTX (
            BLASTN.out.contigs,
            BLASTN.out.nt,
            uniprot_db,
            uniprot_taxids
        )
        ch_versions = ch_versions.mix(DIAMOND_BLASTX.out.versions)
        acdc_contigs = Channel.empty()
        acdc_tax = Channel.empty()

        // BLOBTOOLS
        if (params.blob) {
            BLOBTOOLS (
                DIAMOND_BLASTX.out.contigs,
                DIAMOND_BLASTX.out.nt,
                DIAMOND_BLASTX.out.uniprot,
                DIAMOND_BLASTX.out.real,
                blob_db
            )
            ch_versions = ch_versions.mix(BLOBTOOLS.out.versions)
            acdc_contigs = BLOBTOOLS.out.contigs
            acdc_tax = BLOBTOOLS.out.tax
        }

        if (params.remap) {
            REBLOBTOOLS (
                DIAMOND_BLASTX.out.contigs,
                DIAMOND_BLASTX.out.nt,
                DIAMOND_BLASTX.out.real,
                DIAMOND_BLASTX.out.uniprot,
                blob_db,
                REMAP.out.bam.collect{it[1]},
                REMAP.out.bai.collect{it[1]}
            )
        }

        if (params.acdc) {
            ACDC (
              acdc_contigs,
              acdc_tax,
              kraken_db
            )
        }
    }
    TSNE(ctg)

    // PANGENOME ANALYSIS
    if (params.pangenome) {
        if (params.genusName && params.coreGenesFile) {
            if (params.completeness) {
                COMPLETENESS (
                    ctg,
                    params.genusName,
                    mgpg_db,
                    coreGenesFile
                )
            }
            if (params.tree) {
                TREE (
                    ctg,
                    params.genusName,
                    mgpg_db,
                    coreGenesFile
                )
            }
        }
    }

    faa = Channel.empty()
    prokka_for_split  = Channel.empty()
    ch_multiqc_prokka = Channel.empty()
    if (!euk) {
        PROKKA(ctg, prokka_proteins)
        ch_versions = ch_versions.mix(PROKKA.out.versions)
        PRODIGAL(ctg)
        ch_versions = ch_versions.mix(PRODIGAL.out.versions)
        METARON(ctg, PRODIGAL.out.gff.collect{it[1]})
        ch_versions = ch_versions.mix(METARON.out.versions)
        faa = PROKKA.out.faa
        prokka_for_split  = PROKKA.out.prokka_for_split
        ch_multiqc_prokka = PROKKA.out.prokka_for_split
    } else {
        AUGUSTUS(ctg)
        faa = AUGUSTUS.out.faa
        EUKCC (
            ctg,
            eukcc_db
        )
    }

    if (params.eggnog) {
      EGGNOG (
          faa,
          eggnog_db
      )
      ch_versions = ch_versions.mix(EGGNOG.out.versions)
    }

    // KOFAMSCAN
    kofam_scan = Channel.empty()
    if (params.kofam) {
        KOFAMSCAN (
            faa,
            kofam_profile,
            kofam_kolist
        )
        kofam_scan = KOFAMSCAN.out.txt
        ch_versions = ch_versions.mix(KOFAMSCAN.out.versions)
    }

    // STARAMR
    if (!params.euk) {
        if (params.acquired || params.point) {
            STARAMR (
                ctg,
                params.acquired,
                params.point,
                params.pointfinder_species
            )
            ch_versions = ch_versions.mix(STARAMR.out.versions)
        }
    }

    ch_multiqc_gtdb = Channel.empty()
    if (params.split) {
        split_fa = Channel.empty()
        bin_csv = Channel.empty()
        if (params.split_euk) {
            SPLIT_CHECKM_EUKCC (
                ctg200.collect{it[1]},
                BLOBTOOLS.out.tax_split.collect{it[1]},
                prokka_for_split.collect{it[1]}.ifEmpty([]),
                kofam_scan.collect{it[1]}.ifEmpty([]),
                eukcc_db,
                params.split_bac_level,
                params.split_euk_level
            )
            split_fa = SPLIT_CHECKM_EUKCC.out.fa
            bin_csv = SPLIT_CHECKM_EUKCC.out.csv
        } else {
            SPLIT_CHECKM (
                ctg200.collect{it[1]},
                BLOBTOOLS.out.tax_split.collect{it[1]},
                prokka_for_split.collect{it[1]}.ifEmpty([]),
                kofam_scan.collect{it[1]}.ifEmpty([]),
                params.split_bac_level,
                params.split_euk_level
            )
            split_fa = SPLIT_CHECKM.out.fa
            bin_csv = SPLIT_CHECKM.out.csv
        }

        if (params.graphbin && !params.nanopore) {
            GRAPHBIN (
                contig.collect{it[1]},
                contig_path.collect{it[1]},
                contig_graph.collect{it[1]},
                bin_csv
            )
            ch_versions = ch_versions.mix(GRAPHBIN.out.versions)
        }

        if (params.gtdbtk) {
            GTDBTK (
                split_fa,
                gtdb
            )
            ch_versions = ch_versions.mix(GTDBTK.out.versions)
            ch_multiqc_gtdb = GTDBTK.out.mqc_tsv
        }
    }

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
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_samtools.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_preseq.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_qualimap.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_checkm.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_checkm2.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_gtdb.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_quast.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_prokka.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_kraken.collect{it[1]}.ifEmpty([]))

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
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
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
