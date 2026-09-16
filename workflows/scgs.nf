def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run gongyh/nf-core-scgs --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
    --reads                       Path to input data (must be surrounded with quotes)
    -profile                      Configuration profile to use. Can use multiple (comma separated). Available: conda, docker, singularity, awsbatch, test and more.

    Options:
    --vcf                         Variantion graph construction
    --bulk                        WGS of bulk DNA, not MDA
    --mg                          WGS of bulk DNA, assemble in metagenome mode
    --genome                      Name of iGenomes reference
    --single_end                  Specifies that the input is single end reads
    --snv                         Enable detection of single nucleotide variation
    --cnv                         Enable detection of copy number variation
    --bbmap                       Enable bbmap to remove host-derived contamination
    --doubletd                    Enable detection of doublet
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
    --krona_db                    Krona taxonomy.tab (if used offline)
    --uniprot_db                  Uniprot proteomes database (diamond) !!! time consuming !!!
    --uniprot_taxids              Sequence id to taxa id mapping file
    --kraken2_db                  Kraken2 database
    --kraken1_db                  Kraken1 database (for ACDC)
    --eggnog_db                   EggNOG v4.5.1 database for emapper-1.0.3
    --kofam_profile               KOfam profile database
    --kofam_kolist                KOfam ko_list file
    --augustus_species            Augustus species, default 'saccharomyces'
    --eukcc_db                    EukCC database
    --checkm2_db                  CheckM2 database
    --gtdb                        GTDB database
    --bakta_db                    Bakta database
    --host_ref                    Specify the reference sequence for host removal

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
    --pasa                        Enable PASA scaffolding (default: false)
    --refs_fna                    Genome files for PASA or RAGTAG scaffolding

    Quast options:
    --euk                         Euk genome
    --fungus                      Fungal genome

    Pangenome options:
    --mgpg_db                     Microbiome graph pangenome database
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

def display_header(summary, custom_runName, single_end) {
    // Header log info
    log.info nfcoreHeader()
    //def summary = [:]
    summary['Run Name']         = custom_runName ?: workflow.runName
    summary['Reads']            = params.reads
    summary['Fasta Ref']        = params.fasta
    summary['Data Type']        = single_end ? 'Single-End' : 'Paired-End'
    summary['Bulk']             = params.bulk ? 'Yes' : 'No'
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
// include { BOWTIE2_ALIGN         } from '../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_ALIGN         } from '../modules/local/bowtie2_align'
include { MINIMAP2_ALIGN        } from '../modules/nf-core/minimap2/align/main'
include { QUALIMAP_BAMQC        } from '../modules/nf-core/qualimap/bamqc/main'
include { GENOMAD_ENDTOEND      } from '../modules/nf-core/genomad/endtoend/main'
include { MULTIQC               } from '../modules/nf-core/multiqc/main'

// Import modules from local
include { SAVE_REFERENCE        } from '../modules/local/save_reference'
include { TRIMGALORE            } from '../modules/local/trimgalore'
include { KTUPDATETAXONOMY      } from '../modules/local/ktupdatetaxonomy'
include { KRAKEN                } from '../modules/local/kraken'
include { UMAP                  } from '../modules/local/scanpy/umap'
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
include { SPADES                } from '../modules/local/spades'
include { PANTA; PASA           } from '../modules/local/pasa'
include { COMPLETENESS          } from '../modules/local/pangenome/completeness'
include { TREE                  } from '../modules/local/pangenome/tree'

include { QUAST_REF; QUAST_REF as QUAST_REF0          } from '../modules/local/quast_ref'
include { QUAST_DENOVO; QUAST_DENOVO as QUAST_DENOVO0 } from '../modules/local/quast_denovo'

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
include { BAKTA                 } from '../modules/local/bakta'
include { PRODIGAL              } from '../modules/local/prodigal'
include { UNIOP                 } from '../modules/local/uniop'
include { PROMPREDICT           } from '../modules/local/prompredict'
include { PHISPY                } from '../modules/local/phispy'
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

include { METACOMPASS           } from '../modules/local/metacompass'
include { QUICKMERGE            } from '../modules/local/quickmerge'
include { RAGTAG                } from '../modules/local/ragtag'

/** subworkflow */
include { VG                    } from '../subworkflows/local/vg'

workflow SCGS {
    main:
    /*
 * SET UP CONFIGURATION VARIABLES
 */

// default values
params.reads = "data/*{1,2}.fastq.gz"
params.fasta = false
params.bulk = false
params.outdir = "./results"
params.notrim = false
params.awsregion = "eu-west-1"
params.awsqueue = "default"
params.config_profile_description = null
params.config_profile_contact = null
params.config_profile_url = null
params.email = null
params.maxMultiqcEmailFileSize = 25 * 1024 * 1024
params.genomes = [:]
params.genome = false
params.gff = false
params.vcf = false
params.graph_vcf = null
params.single_end = false
params.fungus = false
params.euk = false
params.genomad_db = null
params.prokka_proteins = null
params.nt_db = null
params.uniprot_db = null
params.uniprot_taxids = null
params.kraken1_db = null
params.kraken2_db = null
params.blob_db = null
params.krona_db = null
params.eggnog_db = null
params.eukcc_db = null
params.checkm2_db = null
params.gtdb = null
params.bakta_db = null
params.mgpg_db = null
params.coreGenesFile = null
params.host_ref = null
params.bbmap = false
params.kofam_profile = null
params.kofam_kolist = null
params.multiqc_config = "$baseDir/assets/multiqc_config.yml"
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0
params.readPaths = null
params.refs_fna = null
params.saveTrimmed = false
params.saveAlignedIntermediates = false
params.kraken = true
params.saturation = false
params.snv = false
params.doubletd = false
params.cnv = false
params.ass = false
params.no_normalize = false
params.mg = false
params.pasa = false
params.genomad = false
params.genus = null
params.checkm2 = true
params.blastn = true
params.evalue = 1e-25
params.blockSize = 2.0
params.blob = true
params.allow_multi_align = false
params.acdc = false
params.pangenome = false
params.genusName = null
params.completeness = false
params.tree = false
params.augustus_species = "saccharomyces"
params.eggnog = true
params.kofam = true
params.acquired = false
params.point = false
params.pointfinder_species = "escherichia_coli"
params.split = false
params.split_euk = false
params.split_bac_level = "genus"
params.split_euk_level = "genus"
params.graphbin = false
params.gtdbtk = true
params.monochrome_logs = false


// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Configurable reference genomes
fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    fasta = file(params.fasta)
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

single_end = params.single_end

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
    def faa = file(params.prokka_proteins)
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
kraken1_db = false
if (params.kraken1_db) {
    kraken1_db = file(params.kraken1_db)
    if( !kraken1_db.exists() ) exit 1, "Kraken v1 database not found: ${params.kraken1_db}"
} else {
    kraken1_db = file("/dev/null")
}

kraken2_db = false
if (params.kraken2_db) {
    kraken2_db = file(params.kraken2_db)
    if( !kraken2_db.exists() ) exit 1, "Kraken v2 database not found: ${params.kraken2_db}"
} else {
    kraken2_db = file("/dev/null")
}

// Configurable Blobtools nodesDB.txt
blob_db = false
if (params.blob_db) {
    blob_db = file(params.blob_db)
    if( !blob_db.exists() ) exit 1, "Blobtools nodesDB.txt not found: ${params.blob_db}"
} else {
    blob_db = file("/dev/null")
}

// Configurable Krona taxonomy.tab
krona_db = false
if (params.krona_db) {
    krona_db = file(params.krona_db)
    if( !krona_db.exists() ) exit 1, "Krona taxonomy.tab not found: ${params.krona_db}"
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

// Configure bakta database
bakta_db = false
if (params.bakta_db) {
    bakta_db  = file(params.bakta_db)
    if ( !bakta_db.exists() ) exit 1, "Bakta database not found: ${params.bakta_db}"
} else {
    bakta_db = file("/dev/null")
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
host_ref = false
if (params.host_ref) {
    host_ref  = file(params.host_ref)
    if ( !host_ref.exists() ) exit 1, "Host reference file not found: ${params.host_ref}"
} else {
    if (params.bbmap) exit 1, "Host reference file not set"
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
ch_multiqc_config = channel.fromPath(params.multiqc_config, checkIfExists: true)
ch_multiqc_custom_config = channel.empty()
ch_multiqc_logo = channel.empty()
ch_output_docs = channel.fromPath("$baseDir/docs/output.md")


// mode
denovo = (params.genome && params.genomes[ params.genome ].bowtie2) || params.fasta ? false : true
bowtie2 = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false


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

if (params.refs_fna) {
    def rfna = file(params.refs_fna, checkIfExists: true)
    if (rfna.size()==1 && rfna.isDirectory()) {
        panta_db = rfna
        refs_fna = channel.empty()
    } else {
        refs_fna = rfna
        panta_db = channel.empty()
    }
} else {
    refs_fna = channel.empty()
    panta_db = channel.empty()
}

summary = [:]

    display_header(summary, custom_runName, single_end)
    ch_versions = channel.empty()

    // FASTQC
    ch_multiqc_fastqc = channel.empty()
    FASTQC ( read_files_fastqc )
    ch_versions       = ch_versions.mix(FASTQC.out.versions)
    ch_multiqc_fastqc = FASTQC.out.zip

    // SAVE_REFERENCE
    if ( params.fasta ) {
        if ( params.gff ) {
            save_reference = SAVE_REFERENCE(fasta, gff)
        } else {
            save_reference = SAVE_REFERENCE(fasta, file("/dev/null"))
        }
    }

    if (bowtie2) {
        bowtie2_index = [bowtie2, file(bowtie2)]
    } else {
        if (params.fasta) {
            def fasta_meta = params.fasta - ~/(\.fasta)?(\.fna)?(\.fa)?$/
            BOWTIE2_BUILD ( [fasta_meta, fasta] )
            bowtie2_index = BOWTIE2_BUILD.out.index
        }
    }
    // TRIM_GALORE
    trimmed_reads = channel.empty()
    ch_multiqc_trim_log = channel.empty()
    ch_multiqc_trim_zip = channel.empty()
    if (params.notrim) {
        if (params.bbmap) {
            bbmap_align = BBMAP_ALIGN(read_files_trimming, host_ref)
            ch_versions = ch_versions.mix(bbmap_align.map { result -> result.versions })
            trimmed_reads = bbmap_align.map { result -> tuple(result.meta, result.clean_fastq) }
        } else {
            trimmed_reads = read_files_trimming
        }
    } else {
        trimgalore = TRIMGALORE(read_files_trimming)
        ch_multiqc_trim_log = trimgalore.map { result -> tuple(result.meta, result.log) }
        ch_multiqc_trim_zip = trimgalore.map { result -> tuple(result.meta, result.zip) }
        ch_versions = ch_versions.mix(trimgalore.map { result -> result.versions })
        if (params.bbmap) {
            bbmap_align = BBMAP_ALIGN(trimgalore.map { result -> tuple(result.meta, result.reads) }, host_ref)
            ch_versions = ch_versions.mix(bbmap_align.map { result -> result.versions })
            trimmed_reads = bbmap_align.map { result -> tuple(result.meta, result.clean_fastq) }
        } else {
            trimmed_reads = trimgalore.map { result -> tuple(result.meta, result.reads) }
        }
    }

    // KRAKEN
    ch_multiqc_kraken = channel.empty()
    if (params.kraken && params.kraken2_db != null) {
        if (!krona_db) {
            krona_download = KTUPDATETAXONOMY()
            krona_db = krona_download.map { result -> result.taxonomy }
        }
        kraken = KRAKEN(
            trimmed_reads,
            kraken2_db,
            krona_db
        )
        UMAP(kraken.map { result -> result.tda }.collect().filter { it -> it.size() >= 4 })
        ch_multiqc_kraken = kraken.map { result -> tuple(result.meta, result.report) }
        ch_versions = ch_versions.mix(kraken.map { result -> result.versions })
    }

    // SATURATION
    if (params.saturation) {
        SATURATION ( trimmed_reads )
    }

    // ALIGN
    if (denovo == false) {
        bowtie2_align = BOWTIE2_ALIGN(
            trimmed_reads,
            bowtie2_index,
            false,
            true
        )
        bb_bam = bowtie2_align.map { result -> tuple(result.meta, result.bam) }
        ch_versions = ch_versions.mix(bowtie2_align.map { result -> result.versions })
    }

    // VG
    if ( params.fasta && params.vcf ) {
        vg = VG (
            fasta,
            trimmed_reads,
            graph_vcf
        )
        ch_versions = ch_versions.mix(vg.ch_versions)
    }

    ch_multiqc_samtools = channel.empty()
    ch_multiqc_preseq   = channel.empty()
    ch_multiqc_qualimap = channel.empty()
    quast_bam = channel.empty()
    quast_bai = channel.empty()
    if ( params.fasta ) {
        ch_samtools_input = bb_bam.combine(save_reference.map { result -> result.bed })
        samtools = SAMTOOLS(ch_samtools_input)
        quast_bam = samtools.map { result -> tuple(result.meta, result.bam) }
        quast_bai = samtools.map { result -> tuple(result.meta, result.bai) }
        ch_samtools_bed = samtools.map { result -> tuple(result.meta, result.bed) }
        ch_versions = ch_versions.mix(samtools.map { result -> result.versions })
        ch_multiqc_samtools = samtools.map { result -> tuple(result.meta, result.stats) }

        preseq = PRESEQ(ch_samtools_bed)
        ch_versions = ch_versions.mix(preseq.map { result -> result.versions })
        ch_multiqc_preseq = preseq.map { result -> tuple(result.meta, result.txt) }

        if ( params.gff ) {
            QUALIMAP_BAMQC (
                quast_bam,
                gff
            )
            ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)
            ch_multiqc_qualimap = QUALIMAP_BAMQC.out.results
        }
        if (params.snv) {
            ch_indelrealign_input = quast_bam.combine(fasta)
            indelrealign = INDELREALIGN(ch_indelrealign_input)
            ch_indelrealign_bam = indelrealign.map { result -> tuple(result.meta, result.bam) }
            ch_indelrealign_bai = indelrealign.map { result -> tuple(result.meta, result.bai) }
            ch_versions = ch_versions.mix(indelrealign.map { result -> result.versions })
        }
        if (!params.bulk && params.snv) {
            monovar = MONOVAR(
                ch_indelrealign_bam.collect { entry -> entry[1] },
                ch_indelrealign_bai.collect { entry -> entry[1] },
                fasta
            )
            ch_versions = ch_versions.mix(monovar.map { result -> result.versions })
            if ( params.doubletd ) {
                doubletd = DOUBLETD(monovar.map { result -> result.vcf })
                ch_versions = ch_versions.mix(doubletd.map { result -> result.versions })
            }
        }
        if (!params.bulk && params.cnv && !single_end) {
            aneufinder = ANEUFINDER(
                quast_bam.collect { entry -> entry[1] },
                quast_bai.collect { entry -> entry[1] }
            )
            ch_versions = ch_versions.mix(aneufinder.map { result -> result.versions })
        }
        ch_circlize_input = ch_samtools_bed.combine(save_reference.map { result -> result.bed })
        circlize = CIRCLIZE(ch_circlize_input)
        ch_versions = ch_versions.mix(circlize.map { result -> result.versions })
    }

    // ASSEMBLY
    ctg200 = channel.empty()
    ctg = channel.empty()
    if ( params.ass ) {
        // NORMALIZE
        if ( params.no_normalize ) {
            trimmed_reads.set{ normalized_reads }
        } else {
            /**
            NORMALIZE(trimmed_reads)
            normalized_reads = NORMALIZE.out.reads
            */
            bbnorm = BBNORM(trimmed_reads)
            normalized_reads = bbnorm.map { result -> tuple(result.meta, result.fastq) }
            ch_versions = ch_versions.mix(bbnorm.map { result -> result.versions })
        }

        spades = SPADES(normalized_reads)
        contig = spades.map { result -> tuple(result.meta, result.contig) }
        contig_path = spades.map { result -> tuple(result.meta, result.contig_path) }
        contig_graph = spades.map { result -> tuple(result.meta, result.contig_graph) }
        ctg200_denovo = spades.map { result -> tuple(result.meta, result.ctg200) }
        ctg_denovo = spades.map { result -> tuple(result.meta, result.ctg) }
        ch_versions = ch_versions.mix(spades.map { result -> result.versions })

        if (params.refs_fna) {
            if (refs_fna.size()>1) {
                panta = PANTA(refs_fna.collect())
                ch_versions = ch_versions.mix(panta.map { result -> result.versions })
                panta_db = panta.map { result -> result.db }
            }
            pasa = PASA(spades.map { result -> tuple(result.meta, result.assembly) }, panta_db)
            ch_versions = ch_versions.mix(pasa.map { result -> result.versions })
            ctg200 = pasa.map { result -> tuple(result.meta, result.ctg200) }
            ctg = pasa.map { result -> tuple(result.meta, result.ctg) }
        } else {
            ctg200 = ctg200_denovo
            ctg = ctg_denovo
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

    // QUAST
    ch_multiqc_quast = channel.empty()
    if (denovo == false) {
        if (params.refs_fna) { // hybrid assembly, add quast for spades
            ch_ctgd_bam_bai = ctg_denovo.join(quast_bam).join(quast_bai).collect(flat: false)
            QUAST_REF0(
                fasta,
                gff,
                ch_ctgd_bam_bai.flatMap { entry -> entry }.map { entry -> entry[1] }.collect(),
                ch_ctgd_bam_bai.flatMap { entry -> entry }.map { entry -> entry[2] }.collect(),
                ch_ctgd_bam_bai.flatMap { entry -> entry }.map { entry -> entry[3] }.collect(),
                euk,
                params.fungus,
                "quast_spades"
            )
        }
        ch_ctg_bam_bai = ctg.join(quast_bam).join(quast_bai).collect(flat: false)
        quast_ref = QUAST_REF(
            fasta,
            gff,
            ch_ctg_bam_bai.flatMap { entry -> entry }.map { entry -> entry[1] }.collect(),
            ch_ctg_bam_bai.flatMap { entry -> entry }.map { entry -> entry[2] }.collect(),
            ch_ctg_bam_bai.flatMap { entry -> entry }.map { entry -> entry[3] }.collect(),
            euk,
            params.fungus,
            "quast_ref"
        )
        ch_multiqc_quast = quast_ref.map { result -> result.tsv }
        ch_versions = ch_versions.mix(quast_ref.map { result -> result.versions })
    } else {
        if (params.refs_fna) { // hybrid assembly, add quast for spades
            QUAST_DENOVO0(
                ctg_denovo.collect { entry -> entry[1] },
                euk,
                params.fungus,
                "quast_spades"
            )
        }
        quast_denovo = QUAST_DENOVO(
            ctg.collect { entry -> entry[1] },
            euk,
            params.fungus,
            "quast_denovo"
        )
        ch_multiqc_quast = quast_denovo.map { result -> result.tsv }
        ch_versions = ch_versions.mix(quast_denovo.map { result -> result.versions })
    }

    // CHECKM_LINEAGEWF
    ch_multiqc_checkm = channel.empty()
    if (!euk) {
        checkm_lineagewf = CHECKM_LINEAGEWF(
            ctg.collect { entry -> entry[1] },
            params.genus ? true : false
        )
        ch_versions = ch_versions.mix(checkm_lineagewf.map { result -> result.versions })
        ch_multiqc_checkm = checkm_lineagewf.map { result -> result.mqc_tsv }
    }

    // CHECKM2
    ch_multiqc_checkm2 = channel.empty()
    if (!euk && params.checkm2 && params.checkm2_db) {
        checkm2 = CHECKM2(
            ctg.collect { entry -> entry[1] },
            'fasta',
            checkm2_db
        )
        ch_versions = ch_versions.mix(checkm2.map { result -> result.versions })
        ch_multiqc_checkm2 = checkm2.map { result -> result.mqc_tsv }
    }

    tax_split = channel.empty()
    if (params.blastn && params.nt_db) {
        // BLASTN
        blastn = BLASTN(
            ctg200,
            nt_db,
            params.evalue
        )
        ch_versions = ch_versions.mix(blastn.map { result -> result.versions })

        // DIAMOND_BLASTS
        diamond_blastx = DIAMOND_BLASTX(
            blastn.map { result -> tuple(result.meta, result.contigs) },
            blastn.map { result -> tuple(result.meta, result.nt) },
            uniprot_db,
            uniprot_taxids,
            params.uniprot_db != null
        )
        ch_versions = ch_versions.mix(diamond_blastx.map { result -> result.versions })
        acdc_contigs = channel.empty()
        acdc_tax = channel.empty()

        // BLOBTOOLS
        if (params.blob && params.blob_db) {
            if (params.no_normalize && !params.refs_fna) {
                ch_blob_input = diamond_blastx.map { result ->
                    tuple(result.meta, result.contigs, result.nt, result.uniprot, result.has_uniprot)
                }
                blobtools = BLOBTOOLS(ch_blob_input, blob_db)
                ch_versions = ch_versions.mix(blobtools.map { result -> result.versions })
                acdc_contigs = blobtools.map { result -> tuple(result.meta, result.contigs) }
                acdc_tax = blobtools.map { result -> tuple(result.meta, result.tax) }
                tax_split = blobtools.map { result -> tuple(result.meta, result.tax_split) }
            } else {
                bowtie2_remap = BOWTIE2_REMAP(ctg200)
                remap_input = trimmed_reads.join(bowtie2_remap.map { result -> tuple(result.meta, result.index) })
                remap = REMAP(remap_input, params.allow_multi_align)
                ch_versions = ch_versions.mix(bowtie2_remap.map { result -> result.versions })
                ch_versions = ch_versions.mix(remap.map { result -> result.versions })
                ch_reblob_input = diamond_blastx
                    .map { result -> tuple(result.meta, result.contigs, result.nt, result.uniprot, result.has_uniprot) }
                    .join(remap.map { result -> tuple(result.meta, result.bam, result.bai) })
                reblobtools = REBLOBTOOLS(ch_reblob_input, blob_db)
                ch_versions = ch_versions.mix(reblobtools.map { result -> result.versions })
                acdc_contigs = reblobtools.map { result -> tuple(result.meta, result.contigs) }
                acdc_tax = reblobtools.map { result -> tuple(result.meta, result.tax) }
                tax_split = reblobtools.map { result -> tuple(result.meta, result.tax_split) }
            }

            if (params.acdc && params.kraken1_db) {
                acdc = ACDC(
                    acdc_contigs,
                    acdc_tax,
                    kraken1_db
                )
                ch_versions = ch_versions.mix(acdc.map { result -> result.versions })
            }
        }
    }
    TSNE(ctg)

    // PANGENOME ANALYSIS
    if (params.pangenome) {
        if (params.genusName && params.coreGenesFile) {
            ch_pangenome_input = ctg.map { meta, contigs ->
                tuple(meta, contigs, params.genusName, mgpg_db, coreGenesFile)
            }
            if (params.completeness) {
                COMPLETENESS(ch_pangenome_input)
            }
            if (params.tree) {
                TREE(ch_pangenome_input)
            }
        }
    }

    faa = channel.empty()
    prokka_for_split  = channel.empty()
    ch_multiqc_prokka = channel.empty()
    if (!euk) {
        prokka = PROKKA(ctg, prokka_proteins)
        ch_versions = ch_versions.mix(prokka.map { result -> result.versions })
        if (params.bakta_db) {
            bakta = BAKTA(
                ctg,
                bakta_db,
                prokka_proteins,
                [] as List<Path>
            )
            ch_versions = ch_versions.mix(bakta.map { result -> result.versions })
        }
        uniop = UNIOP(ctg)
        ch_versions = ch_versions.mix(uniop.map { result -> result.versions })
        prompredict = PROMPREDICT(ctg)
        ch_versions = ch_versions.mix(prompredict.map { result -> result.versions })
        phispy = PHISPY(prokka.map { result -> tuple(result.meta, result.gbk) })
        ch_versions = ch_versions.mix(phispy.map { result -> result.versions })
        faa = prokka.map { result -> tuple(result.meta, result.faa) }
        prokka_for_split = prokka.map { result -> tuple(result.meta, result.prokka_for_split) }
        ch_multiqc_prokka = prokka_for_split
    } else {
        augustus = AUGUSTUS(ctg)
        faa = augustus.map { result -> tuple(result.meta, result.faa) }
        ch_versions = ch_versions.mix(augustus.map { result -> result.versions })
        eukcc = EUKCC(
            ctg,
            eukcc_db
        )
        ch_versions = ch_versions.mix(eukcc.map { result -> result.versions })
    }

    if (params.eggnog && params.eggnog_db) {
        eggnog = EGGNOG(
            faa,
            eggnog_db
        )
        ch_versions = ch_versions.mix(eggnog.map { result -> result.versions })
    }

    // KOFAMSCAN
    kofam_scan = channel.empty()
    if (params.kofam && params.kofam_profile && params.kofam_kolist) {
        kofamscan = KOFAMSCAN(
            faa,
            kofam_profile,
            kofam_kolist
        )
        kofam_scan = kofamscan.map { result -> tuple(result.meta, result.txt) }
        ch_versions = ch_versions.mix(kofamscan.map { result -> result.versions })
    }

    // STARAMR
    if (!params.euk) {
        if (params.acquired || params.point) {
            ch_staramr_input = ctg.map { meta, contigs ->
                tuple(meta, contigs, params.acquired, params.point, params.pointfinder_species ?: '')
            }
            staramr = STARAMR(ch_staramr_input)
            ch_versions = ch_versions.mix(staramr.map { result -> result.versions })
        }
    }

    ch_multiqc_gtdb = channel.empty()
    if (params.split) {
        split_fa = channel.empty()
        bin_csv = channel.empty()
        if (params.split_euk && params.eukcc_db) {
            split_checkm_eukcc = SPLIT_CHECKM_EUKCC(
                ctg200.collect { entry -> entry[1] },
                tax_split.collect { entry -> entry[1] },
                prokka_for_split.collect { entry -> entry[1] }.ifEmpty([]),
                kofam_scan.collect { entry -> entry[1] }.ifEmpty([]),
                eukcc_db,
                params.split_bac_level,
                params.split_euk_level
            )
            split_fa = split_checkm_eukcc.map { result -> result.fa }
            bin_csv = split_checkm_eukcc.map { result -> result.csv }
            ch_versions = ch_versions.mix(split_checkm_eukcc.map { result -> result.versions })
        } else if (!params.split_euk) {
            split_checkm = SPLIT_CHECKM(
                ctg200.collect { entry -> entry[1] },
                tax_split.collect { entry -> entry[1] },
                prokka_for_split.collect { entry -> entry[1] }.ifEmpty([]),
                kofam_scan.collect { entry -> entry[1] }.ifEmpty([]),
                params.split_bac_level,
                params.split_euk_level
            )
            split_fa = split_checkm.map { result -> result.fa }
            bin_csv = split_checkm.map { result -> result.csv }
            ch_versions = ch_versions.mix(split_checkm.map { result -> result.versions })
        }

        if (params.graphbin && !params.refs_fna) {
            graphbin = GRAPHBIN(
                contig.collect { entry -> entry[1] },
                contig_path.collect { entry -> entry[1] },
                contig_graph.collect { entry -> entry[1] },
                bin_csv
            )
            ch_versions = ch_versions.mix(graphbin.map { result -> result.versions })
        }

        if (params.gtdbtk && params.gtdb) {
            gtdbtk = GTDBTK(
                split_fa,
                gtdb
            )
            ch_versions = ch_versions.mix(gtdbtk.map { result -> result.versions })
            ch_multiqc_gtdb = gtdbtk.map { result -> result.mqc_tsv }
        }
    }

    ch_multiqc_versions = channel.empty()
    software_versions = GET_SOFTWARE_VERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    ch_multiqc_versions = software_versions.map { result -> result.mqc_yml }

    // MODULE: MULTIQC
    workflow_summary = create_workflow_summary(summary)
    ch_workflow_summary = channel.value(workflow_summary)

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_fastqc.collect { entry -> entry[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_trim_log.collect { entry -> entry[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_trim_zip.collect { entry -> entry[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_samtools.collect { entry -> entry[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_preseq.collect { entry -> entry[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_qualimap.collect { entry -> entry[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_checkm.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_checkm2.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_gtdb.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_quast.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_prokka.collect { entry -> entry[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_kraken.collect { entry -> entry[1] }.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    OUTPUT_DOCUMENTATION(ch_output_docs)

    emit:
    summary_params = channel.value(summary)
    multiqc_report = MULTIQC.out.report.toList()
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
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  gongyh/nf-core-scgs v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}
