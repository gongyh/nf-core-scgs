def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    nextflow run gongyh/nf-core-scgs --reads '*_R{1,2}.fastq.gz' -profile docker

    Input options:
    --reads <glob>                Input reads glob (default: data/*{1,2}.fastq.gz)
    --readPaths <list>            Structured sample/read list supplied in a Nextflow config
    --single_end                  Treat input reads as single-end
    --bulk                        Process bulk-DNA reads instead of single-cell MDA reads
    --mg                          Enable metagenome assembly mode

    Reference options:
    --genome <name>               Configured iGenomes reference name
    --genomes <map>               Configured iGenomes reference map (normally set in a config)
    --fasta <path>                Reference genome FASTA
    --gff <path>                  Reference genome GFF annotation
    --vcf <path>                  VCF used to construct a variation graph with --fasta
    --euk                         Analyse a eukaryotic genome
    --fungus                      Analyse a fungal genome
    --genus <name>                Supply genus context to CheckM

    Read processing and assembly:
    --notrim                      Skip adapter and quality trimming
    --saveTrimmed                 Publish trimmed reads
    --clip_r1 <int>               Remove bases from the 5' end of read 1
    --clip_r2 <int>               Remove bases from the 5' end of read 2
    --three_prime_clip_r1 <int>   Remove bases from the 3' end of read 1 after trimming
    --three_prime_clip_r2 <int>   Remove bases from the 3' end of read 2 after trimming
    --bbmap                       Remove host-derived reads with BBMap
    --host_ref <path>             Host reference sequence used with --bbmap
    --ass                         Assemble reads with SPAdes
    --no_normalize                Skip read normalization before assembly
    --pasa                        Enable PASA scaffolding
    --refs_fna <path>             Scaffold FASTA files or a PanTA database directory
    --allow_multi_align           Retain secondary and unmapped remapping alignments
    --saveAlignedIntermediates    Publish intermediate alignment BAM files
    --saturation                  Run sequencing saturation analysis

    Variant and genome analyses:
    --snv                         Call single-nucleotide variants with MonoVar
    --cnv                         Call copy-number variants
    --doubletd                    Detect doublets after MonoVar calling
    --genomad                     Run geNomad analysis
    --checkm2                     Run CheckM2 when --checkm2_db is available
    --split                       Split bacterial draft genomes by taxonomic annotation
    --split_euk                   Split eukaryotic draft genomes by taxonomic annotation
    --split_bac_level <level>     Bacterial taxonomic level for --split (default: genus)
    --split_euk_level <level>     Eukaryotic taxonomic level for --split_euk (default: genus)
    --graphbin                    Run GraphBin when splitting de novo assemblies
    --gtdbtk                      Run GTDB-Tk when --gtdb is available
    --pangenome                   Enable pangenome analysis
    --mgpg_db <path>              Microbiome graph pangenome database
    --genusName <name>            Genus name for pangenome analysis
    --coreGenesFile <path>        Core-genes list for pangenome analysis
    --completeness                Calculate pangenome completeness
    --tree                        Build a pangenome phylogenetic tree

    Annotation and ARG analyses:
    --kraken                      Run Kraken2 when --kraken2_db is available
    --blastn                      Run NCBI nt annotation when --nt_db is available
    --blob                        Run BlobTools when --blob_db is available
    --acdc                        Run ACDC when --kraken1_db is available
    --eggnog                      Run EggNOG annotation when --eggnog_db is available
    --kofam                       Run KOfam annotation when profile and KO-list files are available
    --acquired                    Call acquired antimicrobial-resistance genes
    --point                       Call resistance-associated point mutations
    --pointfinder_species <name>  PointFinder species (default: escherichia_coli)
    --evalue <number>             E-value threshold for nt and UniProt searches (default: 1e-25)
    --blockSize <number>          DIAMOND sequence block size in billions of letters (default: 2.0)

    Databases and annotation resources:
    --genomad_db <path>           geNomad database
    --checkm2_db <path>           CheckM2 database
    --nt_db <path>                NCBI nt database
    --blob_db <path>              BlobTools nodesDB.txt
    --krona_db <path>             Krona taxonomy.tab for offline use
    --uniprot_db <path>           UniProt proteomes database
    --uniprot_taxids <path>       UniProt sequence-to-taxonomy mapping
    --kraken2_db <path>           Kraken2 database
    --kraken1_db <path>           Kraken1 database for ACDC
    --eggnog_db <path>            EggNOG database
    --kofam_profile <path>        KOfam profile database
    --kofam_kolist <path>         KOfam KO-list file
    --prokka_proteins <path>      Trusted proteins FASTA for Prokka
    --bakta_db <path>             Bakta database
    --eukcc_db <path>             EukCC database
    --gtdb <path>                 GTDB database
    --augustus_species <name>     Augustus species (default: saccharomyces)

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
    ch_published = channel.empty()

    // FASTQC
    ch_multiqc_fastqc = channel.empty()
    FASTQC ( read_files_fastqc )
    ch_vendor_versions = FASTQC.out.versions
    ch_multiqc_fastqc = FASTQC.out.zip
    ch_published = ch_published.mix(FASTQC.out.html.map { result -> [destination: 'fastqc', files: result] })
    ch_published = ch_published.mix(FASTQC.out.zip.map { result -> [destination: 'fastqc/zips', files: result] })

    // SAVE_REFERENCE
    if ( params.fasta ) {
        if ( params.gff ) {
            save_reference = SAVE_REFERENCE(fasta, gff)
        } else {
            save_reference = SAVE_REFERENCE(fasta, file("/dev/null"))
        }
        if (params.gff) {
            ch_published = ch_published.mix(save_reference.map { result -> [destination: 'reference', files: [result.fa, result.gff, result.bed, result.gc_bed, result.gc_skew_bed, result.genes_bed]] })
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
            trimmed_reads = bbmap_align.map { result -> tuple(result.meta, result.clean_fastq) }
            if (params.saveTrimmed) {
                ch_published = ch_published.mix(bbmap_align.map { result -> [destination: 'remove_hostReads', files: result] })
            }
        } else {
            trimmed_reads = read_files_trimming
        }
    } else {
        trimgalore = TRIMGALORE(read_files_trimming)
        ch_multiqc_trim_log = trimgalore.map { result -> tuple(result.meta, result.logs) }
        ch_multiqc_trim_zip = trimgalore.map { result -> tuple(result.meta, result.fastqc) }
        trimgalore_reads = trimgalore.map { result ->
            def reads = result.meta.single_end ? [result.single_read] : [result.read1, result.read2]
            tuple(result.meta, reads)
        }
        ch_published = ch_published.mix(trimgalore.map { result -> [destination: "trim_galore/${result.meta.id}", files: result.fastqc] })
        ch_published = ch_published.mix(trimgalore.map { result -> [destination: "trim_galore/${result.meta.id}", files: result.logs] })
        if (params.saveTrimmed) {
            ch_published = ch_published.mix(trimgalore_reads.map { _meta, reads -> [destination: 'trim_galore', files: reads] })
        }
        if (params.bbmap) {
            bbmap_align = BBMAP_ALIGN(trimgalore_reads, host_ref)
            trimmed_reads = bbmap_align.map { result -> tuple(result.meta, result.clean_fastq) }
            if (params.saveTrimmed) {
                ch_published = ch_published.mix(bbmap_align.map { result -> [destination: 'remove_hostReads', files: result] })
            }
        } else {
            trimmed_reads = trimgalore_reads
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
        umap = UMAP(kraken.map { result -> result.tda }.collect().filter { it -> it.size() >= 4 })
        ch_multiqc_kraken = kraken.map { result -> tuple(result.meta, result.report) }
        ch_published = ch_published.mix(kraken.map { result -> [destination: 'kraken2', files: result] })
        ch_published = ch_published.mix(umap.map { result -> [destination: 'umap', files: result] })
    }

    // SATURATION
    if (params.saturation) {
        saturation = SATURATION(trimmed_reads)
        ch_published = ch_published.mix(saturation.map { result -> [destination: 'saturation', files: result] })
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
        if (params.saveAlignedIntermediates) {
            ch_published = ch_published.mix(bowtie2_align.map { result -> [destination: 'bowtie2', files: result] })
        }
    }

    // VG
    if ( params.fasta && params.vcf ) {
        vg = VG (
            fasta,
            trimmed_reads,
            graph_vcf
        )
        ch_published = ch_published.mix(vg.published)
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
        ch_multiqc_samtools = samtools.map { result -> tuple(result.meta, result.stats) }
        ch_published = ch_published.mix(samtools.map { result -> [destination: 'bowtie2/stats', files: result.stats] })
        ch_published = ch_published.mix(samtools.map { result -> [destination: 'bowtie2', files: [result.txt, result.pdf]] })
        if (params.saveAlignedIntermediates) {
            ch_published = ch_published.mix(samtools.map { result -> [destination: 'bowtie2', files: [result.bam, result.bai, result.bed, result.versions]] })
        }

        preseq = PRESEQ(ch_samtools_bed)
        ch_multiqc_preseq = preseq.map { result -> tuple(result.meta, result.results) }
        ch_published = ch_published.mix(preseq.map { result -> [destination: '.', files: result.results] })

        if ( params.gff ) {
            QUALIMAP_BAMQC (
                quast_bam,
                gff
            )
            ch_vendor_versions = ch_vendor_versions.mix(QUALIMAP_BAMQC.out.versions)
            ch_multiqc_qualimap = QUALIMAP_BAMQC.out.results
            ch_published = ch_published.mix(QUALIMAP_BAMQC.out.results.map { result -> [destination: 'qualimap_bamqc', files: result] })
            ch_published = ch_published.mix(QUALIMAP_BAMQC.out.versions.map { result -> [destination: 'qualimap_bamqc', files: result] })
        }
        if (params.snv) {
            ch_indelrealign_input = quast_bam.map { meta, bam -> tuple(meta, bam, fasta) }
            indelrealign = INDELREALIGN(ch_indelrealign_input)
            ch_indelrealign_bam = indelrealign.map { result -> tuple(result.meta, result.bam) }
            ch_indelrealign_bai = indelrealign.map { result -> tuple(result.meta, result.bai) }
            ch_published = ch_published.mix(indelrealign.map { result -> [destination: 'gatk', files: result] })
        }
        if (!params.bulk && params.snv) {
            monovar = MONOVAR(
                ch_indelrealign_bam.collect { entry -> entry[1] },
                ch_indelrealign_bai.collect { entry -> entry[1] },
                fasta
            )
            ch_published = ch_published.mix(monovar.map { result -> [destination: 'monovar', files: result.vcf] })
            if ( params.doubletd ) {
                doubletd = DOUBLETD(monovar.map { result -> result.vcf })
                ch_published = ch_published.mix(doubletd.map { result -> [destination: 'doubletd', files: result] })
            }
        }
        if (!params.bulk && params.cnv && !single_end) {
            aneufinder = ANEUFINDER(
                quast_bam.collect { entry -> entry[1] },
                quast_bai.collect { entry -> entry[1] }
            )
            ch_published = ch_published.mix(aneufinder.map { result -> [destination: 'aneufinder', files: result] })
        }
        ch_circlize_input = ch_samtools_bed.combine(save_reference.map { result -> result.bed })
        circlize = CIRCLIZE(ch_circlize_input)
        ch_published = ch_published.mix(circlize.map { result -> [destination: 'circlize', files: result.bed] })
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
            normalized_reads = bbnorm.map { result ->
                def reads = result.meta.single_end ? [result.single_fastq] : [result.fastq1, result.fastq2]
                tuple(result.meta, reads)
            }
        }

        spades = SPADES(normalized_reads)
        ch_published = ch_published.mix(spades.map { result -> [destination: 'spades', files: result] })
        contig = spades.map { result -> tuple(result.meta, result.contig) }
        contig_path = spades.map { result -> tuple(result.meta, result.contig_path) }
        contig_graph = spades.map { result -> tuple(result.meta, result.contig_graph) }
        ctg200_denovo = spades.map { result -> tuple(result.meta, result.ctg200) }
        ctg_denovo = spades.map { result -> tuple(result.meta, result.ctg) }

        if (params.refs_fna) {
            if (refs_fna.size()>1) {
                panta = PANTA(refs_fna.collect())
                panta_db = panta.map { result -> result.db }
                ch_published = ch_published.mix(panta.map { result -> [destination: 'pasa', files: result] })
            }
            pasa = PASA(spades.map { result -> tuple(result.meta, result.assembly) }, panta_db)
            ch_published = ch_published.mix(pasa.map { result -> [destination: 'pasa', files: result] })
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
        ch_vendor_versions = ch_vendor_versions.mix(GENOMAD_ENDTOEND.out.versions)
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.aggregated_classification.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.taxonomy.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.provirus.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.compositions.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.calibrated_classification.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.plasmid_fasta.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.plasmid_genes.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.plasmid_proteins.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.plasmid_summary.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.virus_fasta.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.virus_genes.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.virus_proteins.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.virus_summary.map { result -> [destination: 'genomad', files: result] })
        ch_published = ch_published.mix(GENOMAD_ENDTOEND.out.versions.map { result -> [destination: 'genomad', files: result] })
    }

    // QUAST
    ch_multiqc_quast = channel.empty()
    if (denovo == false) {
        if (params.refs_fna) { // hybrid assembly, add quast for spades
            ch_ctgd_bam_bai = ctg_denovo.join(quast_bam).join(quast_bai).collect(flat: false)
            quast_ref0 = QUAST_REF0(
                fasta,
                gff,
                ch_ctgd_bam_bai.flatMap { entry -> entry }.map { entry -> entry[1] }.collect(),
                ch_ctgd_bam_bai.flatMap { entry -> entry }.map { entry -> entry[2] }.collect(),
                ch_ctgd_bam_bai.flatMap { entry -> entry }.map { entry -> entry[3] }.collect(),
                euk,
                params.fungus,
                "quast_spades"
            )
            ch_published = ch_published.mix(quast_ref0.map { result -> [destination: 'quast', files: result] })
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
        ch_published = ch_published.mix(quast_ref.map { result -> [destination: 'quast', files: result] })
    } else {
        if (params.refs_fna) { // hybrid assembly, add quast for spades
            quast_denovo0 = QUAST_DENOVO0(
                ctg_denovo.collect { entry -> entry[1] },
                euk,
                params.fungus,
                "quast_spades"
            )
            ch_published = ch_published.mix(quast_denovo0.map { result -> [destination: 'quast', files: result] })
        }
        quast_denovo = QUAST_DENOVO(
            ctg.collect { entry -> entry[1] },
            euk,
            params.fungus,
            "quast_denovo"
        )
        ch_multiqc_quast = quast_denovo.map { result -> result.tsv }
        ch_published = ch_published.mix(quast_denovo.map { result -> [destination: 'quast', files: result] })
    }

    // CHECKM_LINEAGEWF
    ch_multiqc_checkm = channel.empty()
    if (!euk) {
        checkm_lineagewf = CHECKM_LINEAGEWF(
            ctg.collect { entry -> entry[1] },
            params.genus ? true : false
        )
        ch_multiqc_checkm = checkm_lineagewf.map { result -> result.mqc_tsv }
        ch_published = ch_published.mix(checkm_lineagewf.map { result -> [destination: 'CheckM', files: result] })
    }

    // CHECKM2
    ch_multiqc_checkm2 = channel.empty()
    if (!euk && params.checkm2 && params.checkm2_db) {
        checkm2 = CHECKM2(
            ctg.collect { entry -> entry[1] },
            'fasta',
            checkm2_db
        )
        ch_multiqc_checkm2 = checkm2.map { result -> result.mqc_tsv }
        ch_published = ch_published.mix(checkm2.map { result -> [destination: 'CheckM2', files: result] })
    }

    tax_split = channel.empty()
    if (params.blastn && params.nt_db) {
        // BLASTN
        blastn = BLASTN(
            ctg200,
            nt_db,
            Float.valueOf(params.evalue.toString())
        )
        ch_published = ch_published.mix(blastn.map { result -> [destination: 'blob', files: result] })

        // DIAMOND_BLASTS
        diamond_blastx = DIAMOND_BLASTX(
            blastn.map { result -> tuple(result.meta, result.contigs) },
            blastn.map { result -> tuple(result.meta, result.nt) },
            uniprot_db,
            uniprot_taxids,
            params.uniprot_db != null
        )
        ch_published = ch_published.mix(diamond_blastx.map { result -> [destination: 'blob', files: result] })
        acdc_contigs = channel.empty()
        acdc_tax = channel.empty()

        // BLOBTOOLS
        if (params.blob && params.blob_db) {
            if (params.no_normalize && !params.refs_fna) {
                ch_blob_input = diamond_blastx.map { result ->
                    tuple(result.meta, result.contigs, result.nt, result.uniprot, result.has_uniprot)
                }
                blobtools = BLOBTOOLS(ch_blob_input, blob_db)
                acdc_contigs = blobtools.map { result -> tuple(result.meta, result.contigs) }
                acdc_tax = blobtools.map { result -> tuple(result.meta, result.tax) }
                tax_split = blobtools.map { result -> tuple(result.meta, result.tax_split) }
                ch_published = ch_published.mix(blobtools.map { result -> [destination: 'blob', files: result] })
            } else {
                bowtie2_remap = BOWTIE2_REMAP(ctg200)
                remap_input = trimmed_reads.join(bowtie2_remap.map { result -> tuple(result.meta, result.index) })
                remap = REMAP(remap_input, params.allow_multi_align)
                ch_published = ch_published.mix(remap.map { result -> [destination: 'remap', files: result] })
                ch_reblob_input = diamond_blastx
                    .map { result -> tuple(result.meta, result.contigs, result.nt, result.uniprot, result.has_uniprot) }
                    .join(remap.map { result -> tuple(result.meta, result.bam, result.bai) })
                reblobtools = REBLOBTOOLS(ch_reblob_input, blob_db)
                acdc_contigs = reblobtools.map { result -> tuple(result.meta, result.contigs) }
                acdc_tax = reblobtools.map { result -> tuple(result.meta, result.tax) }
                tax_split = reblobtools.map { result -> tuple(result.meta, result.tax_split) }
                ch_published = ch_published.mix(reblobtools.map { result -> [destination: 'reblob', files: result] })
            }

            if (params.acdc && params.kraken1_db) {
                acdc = ACDC(
                    acdc_contigs,
                    acdc_tax,
                    kraken1_db
                )
                ch_published = ch_published.mix(acdc.map { result -> [destination: 'acdc', files: result] })
            }
        }
    }
    tsne = TSNE(ctg)
    ch_published = ch_published.mix(tsne.map { result -> [destination: 'tsne', files: result] })

    // PANGENOME ANALYSIS
    if (params.pangenome) {
        if (params.genusName && params.coreGenesFile) {
            ch_pangenome_input = ctg.map { meta, contigs ->
                tuple(meta, contigs, params.genusName, mgpg_db, coreGenesFile)
            }
            if (params.completeness) {
                completeness = COMPLETENESS(ch_pangenome_input)
                ch_published = ch_published.mix(completeness.map { result -> [destination: 'mgpg', files: result] })
            }
            if (params.tree) {
                tree = TREE(ch_pangenome_input)
                ch_published = ch_published.mix(tree.map { result -> [destination: 'mgpg', files: result] })
            }
        }
    }

    faa = channel.empty()
    prokka_for_split  = channel.empty()
    ch_multiqc_prokka = channel.empty()
    if (!euk) {
        prokka = PROKKA(ctg, prokka_proteins)
        ch_published = ch_published.mix(prokka.map { result -> [destination: 'prokka', files: result] })
        if (params.bakta_db) {
            bakta = BAKTA(
                ctg,
                bakta_db,
                prokka_proteins,
                [] as List<Path>
            )
            ch_published = ch_published.mix(bakta.map { result -> [destination: 'bakta', files: result] })
        }
        uniop = UNIOP(ctg)
        prompredict = PROMPREDICT(ctg)
        phispy = PHISPY(prokka.map { result -> tuple(result.meta, result.gbk) })
        ch_published = ch_published.mix(uniop.map { result -> [destination: 'operons', files: result] })
        ch_published = ch_published.mix(prompredict.map { result -> [destination: 'promoters', files: result] })
        ch_published = ch_published.mix(phispy.map { result -> [destination: 'prophages', files: result] })
        faa = prokka.map { result -> tuple(result.meta, result.faa) }
        prokka_for_split = prokka.map { result -> tuple(result.meta, result.prokka_for_split) }
        ch_multiqc_prokka = prokka_for_split
    } else {
        augustus = AUGUSTUS(ctg)
        ch_published = ch_published.mix(augustus.map { result -> [destination: 'augustus', files: result] })
        faa = augustus.map { result -> tuple(result.meta, result.faa) }
        eukcc = EUKCC(
            ctg,
            eukcc_db
        )
        ch_published = ch_published.mix(eukcc.map { result -> [destination: 'eukcc', files: result] })
    }

    if (params.eggnog && params.eggnog_db) {
        eggnog = EGGNOG(
            faa,
            eggnog_db
        )
        ch_published = ch_published.mix(eggnog.map { result -> [destination: 'eggnog', files: result] })
    }

    // KOFAMSCAN
    kofam_scan = channel.empty()
    if (params.kofam && params.kofam_profile && params.kofam_kolist) {
        kofamscan = KOFAMSCAN(
            faa,
            kofam_profile,
            kofam_kolist
        )
        kofam_scan = kofamscan.flatMap { result -> result.txt }
        ch_published = ch_published.mix(kofamscan.map { result -> [destination: 'kofam', files: result] })
    }

    // STARAMR
    if (!params.euk) {
        if (params.acquired || params.point) {
            ch_staramr_input = ctg.map { meta, contigs ->
                tuple(meta, contigs, params.acquired, params.point, params.pointfinder_species ?: '')
            }
            staramr = STARAMR(ch_staramr_input)
            ch_published = ch_published.mix(staramr.map { result -> [destination: 'ARG', files: result] })
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
                kofam_scan.collect().ifEmpty([]),
                eukcc_db,
                params.split_bac_level,
                params.split_euk_level
            )
            split_fa = split_checkm_eukcc.map { result -> result.fa }
            bin_csv = split_checkm_eukcc.map { result -> result.csv }
            ch_published = ch_published.mix(split_checkm_eukcc.map { result -> [destination: '.', files: result] })
        } else if (!params.split_euk) {
            split_checkm = SPLIT_CHECKM(
                ctg200.collect { entry -> entry[1] },
                tax_split.collect { entry -> entry[1] },
                prokka_for_split.collect { entry -> entry[1] }.ifEmpty([]),
                kofam_scan.collect().ifEmpty([]),
                params.split_bac_level,
                params.split_euk_level
            )
            split_fa = split_checkm.map { result -> result.fa }
            bin_csv = split_checkm.map { result -> result.csv }
            ch_published = ch_published.mix(split_checkm.map { result -> [destination: '.', files: result] })
        }

        if (params.graphbin && !params.refs_fna) {
            graphbin = GRAPHBIN(
                contig.collect { entry -> entry[1] },
                contig_path.collect { entry -> entry[1] },
                contig_graph.collect { entry -> entry[1] },
                bin_csv
            )
            ch_published = ch_published.mix(graphbin.map { result -> [destination: 'graphbin', files: result] })
        }

        if (params.gtdbtk && params.gtdb) {
            gtdbtk = GTDBTK(
                split_fa,
                gtdb
            )
            ch_multiqc_gtdb = gtdbtk.map { result -> result.mqc_tsv }
            ch_published = ch_published.mix(gtdbtk.map { result -> [destination: 'gtdb', files: result] })
        }
    }

    ch_multiqc_versions = channel.empty()
    software_versions = GET_SOFTWARE_VERSIONS(
        channel.topic('local_versions')
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
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  gongyh/nf-core-scgs v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}
