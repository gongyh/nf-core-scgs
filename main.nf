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


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
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
      --bulk                        WGS of bulk DNA, not MDA
      --nanopore                    Nanopore sequencing data, force single_end, assembly using Canu
      --genome                      Name of iGenomes reference
      --single_end                  Specifies that the input is single end reads
      --snv                         Enable detection of single nucleotide variation
      --cnv                         Enable detection of copy number variation
      --saturation                  Enable sequencing saturation analysis
      --ass                         Assemble using SPAdes
      --split                       Split the draft genomes and annotation

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
      --three_prime_clip_r2 [int]   Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed

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
params.notrim = false
params.saveTrimmed = false
params.allow_multi_align = false
params.saveAlignedIntermediates = false
params.no_normalize = false
params.euk = false
params.fungus = false
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

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// TODO nf-core: Add any reference files that are needed
// Configurable reference genomes
fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if ( params.fasta ){
    fasta = file(params.fasta)
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

/*
 * Create a channel for input read files
 */
if(params.readPaths){
    if(single_end){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_files_fastqc; read_files_trimming }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_files_fastqc; read_files_trimming }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: single_end ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .into { read_files_fastqc; read_files_trimming }
}


// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
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


/*
 * Parse software version numbers
 */
process get_software_versions {
    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml1, software_versions_yaml2

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    trim_galore --version &> v_trim_galore.txt
    bowtie2 --version &> v_bowtie2.txt
    minimap2 -V &> v_minimap2.txt
    samtools --version &> v_samtools.txt
    bedtools --version &> v_bedtools.txt
    preseq &> v_preseq.txt
    qualimap -h &> v_qualimap.txt
    if [ x=`picard MarkDuplicates --version &> v_picard.txt` = 1 ]; then return 0; fi
    #gatk3 -version &> v_gatk.txt
    Rscript -e 'print(packageVersion("AneuFinder"))' &> v_AneuFinder.txt
    spades.py --version &> v_spades.txt
    canu -version &> v_canu.txt
    blastn -version &> v_blast.txt
    quast.py --version &> v_quast.txt
    multiqc --version &> v_multiqc.txt
    diamond version &> v_diamond.txt
    kraken --version | grep Kraken &> v_kraken.txt
    head -n 1 /opt/conda/envs/nf-core-gongyh-scgs/lib/python3.6/site-packages/checkm/VERSION &> v_checkm.txt
    prokka -v &> v_prokka.txt
    emapper.py --version | grep emapper &> v_eggnogmapper.txt
    echo 'v0.0.1' > v_monovar.txt
    blobtools -v &> v_blobtools.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * Store reference
 */
process save_reference {
    publishDir path: "${params.outdir}/reference_genome", mode: 'copy'

    input:
    file fasta from fasta
    file gff from gff

    output:
    file "genome.fa"
    file "genome.gff"
    file "*.bed"
    file "genome.bed" into genome_circlize, genome_samtools

    when:
    params.fasta && params.gff

    script:
    """
    ln -s ${fasta} genome.fa
    ln -s ${gff} genome.gff
    fa2bed.py genome.fa
    cat genome.gff | grep \$'\tgene\t' | bedtools sort | cut -f1,4,5,7 > genes.bed
    """
}

/*
 * STEP 0 - Split Build Bowtie2 index if necessary
 */

bowtie2 = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false
if(bowtie2){
    bowtie2_index = file(bowtie2)
} else {
process prepare_bowtie2 {
    publishDir path: "${params.outdir}/reference_genome", mode: 'copy'

    input:
    file fasta from fasta

    output:
    file "Bowtie2Index" into bowtie2_index

    when:
    params.fasta

    script:
    """
    mkdir -p Bowtie2Index; cd Bowtie2Index
    ln -s ../${fasta} genome.fa
    bowtie2-build genome.fa genome
    """
}
}

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results1, fastqc_results2

    script:
    """
    fastqc --threads ${task.cpus} -q $reads
    """
}

/*
 * STEP 2 - Trim Galore!
 */
if(params.notrim){
    read_files_trimming.map {name, reads -> reads}
        .into { trimmed_reads; trimmed_reads_for_spades; trimmed_reads_for_kraken; trimmed_reads_for_kmer }
    trimgalore_results1 = file('/dev/null')
    trimgalore_results2 = file('/dev/null')
    trimgalore_fastqc_reports1 = file('/dev/null')
    trimgalore_fastqc_reports2 = file('/dev/null')
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(name), file(reads) from read_files_trimming

        output:
        file '*.fq.gz' into trimmed_reads, trimmed_reads_for_spades, trimmed_reads_for_kraken, trimmed_reads_for_kmer
        file '*trimming_report.txt' into trimgalore_results1, trimgalore_results2
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports1, trimgalore_fastqc_reports2

        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        if (single_end) {
            """
            trim_galore --trim-n --max_n 0 --fastqc --gzip --fastqc_args \"--threads ${task.cpus}\" --cores 4 $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --trim-n --max_n 0 --fastqc --gzip --fastqc_args \"--threads ${task.cpus}\" --cores 4 $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}

/*
 * STEP 2.1 - kraken
 */
process kraken {
    tag "$prefix"
    publishDir path: "${params.outdir}/kraken", mode: 'copy'

    input:
    file reads from trimmed_reads_for_kraken
    file db from kraken_db

    output:
    file "${prefix}.report" into kraken_for_mqc1, kraken_for_mqc2
    file "${prefix}.krona.html"

    when:
    params.kraken_db

    script:
    prefix = reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    mode = single_end ? "" : "--paired"
    """
    kraken -db $db --threads ${task.cpus} --fastq-input --gzip-compressed ${mode} --check-names --output ${prefix}.krk $reads
    kraken-report -db $db ${prefix}.krk > ${prefix}.report
    cut -f2,3 ${prefix}.krk > ${prefix}.f23
    ktImportTaxonomy -o ${prefix}.krona.html ${prefix}.f23
    """
}

/*
 * STEP 2.2 - saturation
 */
process saturation {
    tag "$prefix"
    publishDir path: "${params.outdir}/saturation", mode: 'copy'

    input:
    file reads from trimmed_reads_for_kmer

    output:
    file "${prefix}_kmer.pdf"
    file "${prefix}_cov31_*.csv"

    when:
    params.saturation

    script:
    prefix = reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    R1 = reads[0].toString()
    if (single_end) {
    """
    fastp -i $R1 -A -G -Q -L -s 10 -d 0 -o ${prefix}_split.fq.gz
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
    R2 = reads[1].toString()
    """
    fastp -i $R1 -I $R2 -A -G -Q -L -s 10 -d 0 -o ${prefix}_split_R1.fq.gz -O ${prefix}_split_R2.fq.gz
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

/*
 * STEP 3 - align with bowtie2
 */
if ( params.nanopore ) {

process minimap2 {
    tag "$prefix"
    publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/bowtie2" : params.outdir }, mode: 'copy',
               saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }
    input:
    file reads from trimmed_reads
    file fasta from fasta
    
    output:
    file '*.bam' into bb_bam
    
    when:
    denovo == false
    
    script:
    prefix = reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    filtering = params.allow_multi_align ? '' : "-F 256"
    R1 = reads[0].toString()
    """
    minimap2 -x map-ont -a $fasta $R1 | samtools view -b -q 40 -F 4 $filtering - > ${prefix}.bam
    """
}

} else {

process bowtie2 {
    tag "$prefix"
    publishDir path: { params.saveAlignedIntermediates ? "${params.outdir}/bowtie2" : params.outdir }, mode: 'copy',
               saveAs: {filename -> params.saveAlignedIntermediates ? filename : null }

    input:
    file reads from trimmed_reads
    file index from bowtie2_index

    output:
    file '*.bam' into bb_bam

    when:
    denovo == false

    script:
    prefix = reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    filtering = params.allow_multi_align ? '' : "| samtools view -b -q 40 -F 4 -F 256 -"
    R1 = reads[0].toString()
    if (single_end) {
    """
    bowtie2 -x ${index}/genome -p ${task.cpus} -U $R1 | samtools view -bT $index - $filtering > ${prefix}.bam
    """
    } else {
    R2 = reads[1].toString()
    """
    bowtie2 --no-mixed --no-discordant -X 1000 -x ${index}/genome -p ${task.cpus} -1 $R1 -2 $R2 | samtools view -bT $index - $filtering > ${prefix}.bam
    """
    }
}

}

/*
 * STEP 4 - post-alignment processing
 */
process samtools {
    tag "${prefix}"
    publishDir path: "${pp_outdir}", mode: 'copy',
               saveAs: { filename ->
                   if (filename.indexOf(".stats.txt") > 0) "stats/$filename"
                   else if (filename.indexOf("_bins.txt") > 0) filename
                   else if (filename.indexOf("_pdrc.pdf") > 0) filename
                   else params.saveAlignedIntermediates ? filename : null
               }

    input:
    file bam from bb_bam
    file genome from genome_samtools

    output:
    file '*.markdup.bam' into bam_for_qualimap, bam_for_aneufinder, bam_for_realign, bam_for_quast
    file '*.markdup.bam.bai' into bai_for_qualimap, bai_for_aneufinder, bai_for_realign, bai_for_quast
    file '*.markdup.bed' into bed_for_circlize, bed_for_preseq
    file '*.stats.txt' into samtools_stats
    file "${prefix}_1k_bins.txt"
    file "${prefix}_pdrc.pdf"

    script:
    pp_outdir = "${params.outdir}/bowtie2"
    prefix = bam.baseName
    """
    samtools sort -o ${prefix}.sorted.bam $bam
    samtools index ${prefix}.sorted.bam
    picard MarkDuplicates I=${prefix}.sorted.bam O=${prefix}.markdup.bam M=metrics.txt AS=true
    samtools index ${prefix}.markdup.bam
    bedtools bamtobed -i ${prefix}.markdup.bam | sort -T /tmp -k 1,1 -k 2,2n -k 3,3n -k 6,6 > ${prefix}.markdup.bed
    samtools stats -t ${genome} ${prefix}.markdup.bam > ${prefix}.stats.txt
    # uniformity
    cut -f1,3 ${genome} > ref.genome
    genomeCoverageBed -ibam ${prefix}.markdup.bam -d -g ref.genome > ${prefix}.raw.cov
    meanCov=\$(awk 'BEGIN{ total=0; base=0 } { total=total+\$3; base=base+1 } END{ printf total/base }' ${prefix}.raw.cov)
    awk -v mc=\$meanCov -F'\t' '{print \$1"\t"\$2"\t"\$3/mc}' ${prefix}.raw.cov > ${prefix}.relative.cov
    awk '{sum+=\$3} (NR%1000)==0{print sum/1000; sum=0;}' ${prefix}.relative.cov > ${prefix}_1k_bins.txt
    plotProp.R ${prefix}_1k_bins.txt ${prefix}
    """
}

/*
 * STEP 4.1 - predicting library complexity and genome coverage using preseq
 */
process preseq {
    tag "${prefix}"
    publishDir path: "${pp_outdir}", mode: 'copy',
               saveAs: { filename ->
                   if (filename.indexOf(".txt") > 0) filename
                   else if (filename.indexOf(".pdf") > 0) filename
                   else null }

    input:
    file sbed from bed_for_preseq

    output:
    file '*.txt' into preseq_for_multiqc
    file '*.pdf'

    when:
    !params.nanopore

    script:
    pp_outdir = "${params.outdir}/preseq"
    prefix = sbed.toString() - ~/(\.markdup\.bed)?(\.markdup)?(\.bed)?$/
    mode = single_end ? "" : "-P"
    if (params.bulk) {
    """
    preseq c_curve ${mode} -s 1e+5 -o ${prefix}_c.txt $sbed
    preseq lc_extrap ${mode} -s 1e+5 -D -o ${prefix}_lc.txt $sbed
    plotPreSeq.R ${prefix}_lc.txt ${prefix}_lc
    """
    } else {
    """
    preseq c_curve ${mode} -s 1e+5 -o ${prefix}_c.txt $sbed
    preseq lc_extrap ${mode} -s 1e+5 -D -o ${prefix}_lc.txt $sbed
    plotPreSeq.R ${prefix}_lc.txt ${prefix}_lc
    preseq gc_extrap -w 1000 -s 1e+7 -D -o ${prefix}_gc.txt $sbed
    plotPreSeq.R ${prefix}_gc.txt ${prefix}_gc
    """
    }
}

/*
 * STEP 4.2 - quality control of alignment sequencing data using QualiMap
 */
process qualimap {
    publishDir path: "${pp_outdir}", mode: 'copy'

    input:
    file ("*") from bam_for_qualimap.collect()
    file ("*") from bai_for_qualimap.collect()
    file gff from gff

    output:
    file '*.markdup_stats' into qualimap_for_multiqc
    file 'multi-bamqc'

    script:
    pp_outdir = "${params.outdir}/qualimap"
    """
    ls *.markdup.bam > bams.txt
    let num=`ls *.bam | wc -l`
    if [ \$num == 1 ]; then
      qualimap bamqc -c -bam *.markdup.bam -gff $gff -outdir multi-bamqc
      ln -s multi-bamqc Sample.markdup_stats
    else
      cat bams.txt | awk '{split(\$1,a,".markdup.bam"); print a[1]"\t"\$1}' > inputs.txt
      JAVA_MEM_SIZE=${task.memory.toGiga()}G qualimap multi-bamqc -r -c -d inputs.txt -gff $gff -outdir multi-bamqc
    fi
    """
}

/*
 * STEP 4.3.0 - Realign InDels
 */
process IndelRealign {
    tag "${prefix}"
    publishDir path: "${pp_outdir}", mode: 'copy'

    input:
    file bam from bam_for_realign
    file fa from fasta

    output:
    file '*.realign.bam' into bam_for_monovar
    file '*.realign.bam.bai' into bai_for_monovar

    when:
    params.snv && !params.nanopore

    script:
    pp_outdir = "${params.outdir}/gatk"
    prefix = bam.toString() - ~/(\.markdup\.bam)?(\.markdup)?(\.bam)?$/
    """
    samtools faidx $fa
    picard CreateSequenceDictionary R=$fa
    picard AddOrReplaceReadGroups I=$bam O=${prefix}_rg.bam RGLB=lib RGPL=illumina RGPU=run RGSM=${prefix}
    samtools index ${prefix}_rg.bam
    gatk3 -T RealignerTargetCreator -R $fa -I ${prefix}_rg.bam -o indels.intervals
    gatk3 -T IndelRealigner -R $fa -I ${prefix}_rg.bam -targetIntervals indels.intervals -o ${prefix}.realign.bam
    #java -Xmx4g -jar ${workflow.projectDir}/bin/srma-0.1.15.jar I=${prefix}_rg.bam O=${prefix}.realign.bam R=${fa}
    samtools index ${prefix}.realign.bam
    """
}

/*
 * STEP 4.3 - SNV detection using MonoVar
 */
process monovar {
    publishDir path: "${pp_outdir}", mode: 'copy',
               saveAs: { filename ->
                   if (filename.indexOf(".vcf") > 0) "$filename" else null }

    input:
    file ("*") from bam_for_monovar.collect()
    file ("*") from bai_for_monovar.collect()
    file fa from fasta

    output:
    file 'monovar.vcf' into monovar_vcf

    when:
    !params.bulk && params.snv && !params.nanopore

    script:
    pp_outdir = "${params.outdir}/monovar"
    """
    ls *.bam > bams.txt
    samtools mpileup -B -d 10000 -q 40 -f $fa -b bams.txt | /opt/MonoVar/src/monovar.py -f $fa -o monovar.vcf -m ${task.cpus} -b bams.txt
    """
}

/*
 * STEP 4.4 - CNV detection using AneuFinder
 */
process aneufinder {
    publishDir path: "${pp_outdir}", mode: 'copy'

    input:
    file ("bams/*") from bam_for_aneufinder.collect()
    file ("bams/*") from bai_for_aneufinder.collect()

    output:
    file 'CNV_output' into cnv_output

    when:
    !params.bulk && params.cnv && !single_end && !params.nanopore

    script:
    pp_outdir = "${params.outdir}/aneufinder"
    """
    aneuf.R ./bams CNV_output ${task.cpus}
    """
}

/*
 * STEP 5 - Prepare files for Circlize
 */
process circlize {
    tag "${prefix}"
    publishDir "${params.outdir}/circlize", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".bed") > 0) "$filename" else null
            }

    input:
    file sbed from bed_for_circlize
    file refbed from genome_circlize

    output:
    file "${prefix}-cov200.bed"

    shell:
    prefix = sbed.toString() - ~/(\.markdup\.bed)?(\.markdup)?(\.bed)?$/
    """
    bedtools makewindows -b $refbed -w 200 > genome.200.bed
    bedtools coverage -mean -b $sbed -a genome.200.bed | sort -k 1V,1 -k 2n,2 -k 3n,3 > ${prefix}-cov200.bed
    """
}

/*
 * STEP 6.0 - Normalize reads for assembly
 */
if (params.no_normalize) {
    trimmed_reads_for_spades.set { normalized_reads_for_assembly }
} else {

process normalize {
    tag "${prefix}"

    input:
    file clean_reads from trimmed_reads_for_spades

    output:
    file "*_norm*.fastq.gz" into normalized_reads_for_assembly

    when:
    params.ass

    script:
    prefix = clean_reads[0].toString() - ~/(\.R1)?(_1)?(_R1)?(_trimmed)?(_combined)?(\.1_val_1)?(_1_val_1)?(_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?(\.bz2)?$/
    R1 = clean_reads[0].toString()
    mode = params.bulk ? "bulk" : "mda"
    if (single_end) {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
      ln -s $R1 ${prefix}_norm.fastq.gz
    else
      normalize-by-median.py -k 31 -C 40 --gzip -M 4e+9 -R ${prefix}_norm.report -o ${prefix}_norm.fastq.gz $R1
    fi
    """
    } else {
    R2 = clean_reads[1].toString()
    """
    if [ \"${mode}\" == \"bulk\" ]; then
      ln -s $R1 ${prefix}_norm_R1.fastq.gz
      ln -s $R2 ${prefix}_norm_R2.fastq.gz
    else
      gzip -cd $R1 | fastx_renamer -n COUNT -i /dev/stdin -Q33 -z -o ${prefix}_rename_R1_fq.gz
      gzip -cd $R2 | fastx_renamer -n COUNT -i /dev/stdin -Q33 -z -o ${prefix}_rename_R2_fq.gz
      interleave-reads.py ${prefix}_rename_R1_fq.gz ${prefix}_rename_R2_fq.gz | normalize-by-median.py -k 31 -C 40 -M 4e+9 -p --gzip -R ${prefix}_norm.report -o ${prefix}_nbm.fastq.gz /dev/stdin
      split-paired-reads.py -1 ${prefix}_norm_R1.fastq.gz -2 ${prefix}_norm_R2.fastq.gz --gzip ${prefix}_nbm.fastq.gz
    fi
    """
    }
}

}


/*
 * STEP 6 - Assemble using SPAdes
 */
if ( params.nanopore ) {

process canu {
    tag "${prefix}"
    publishDir path: "${params.outdir}/spades", mode: 'copy'

    input:
    file clean_reads from normalized_reads_for_assembly

    output:
    file "${prefix}.ctg200.fasta" into contigs_for_nt, contigs_for_split
    file "${prefix}.ctgs.fasta" into contigs_for_quast1, contigs_for_quast2, contigs_for_checkm, contigs_for_prokka, contigs_for_prodigal, contigs_for_resfinder, contigs_for_pointfinder, contigs_for_tsne, contigs_for_augustus, contigs_for_eukcc

    when:
    params.ass

    script:
    prefix = clean_reads[0].toString() - ~/(_trimmed)?(_norm)?(_combined)?(\.R1)?(_1)?(_R1)?(\.1_val_1)?(_1_val_1)?(_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?(\.bz2)?$/
    R1 = clean_reads[0].toString()
    mode = params.bulk ? "bulk" : "mda"
    """
    if [ \"${mode}\" == \"bulk\" ]; then
      #canu -d ${prefix}.spades_out -p ${prefix} genomeSize=4m useGrid=false maxThreads=${task.cpus} maxMemory=${task.memory.toGiga()}g -nanopore $R1
      flye --nano-raw $R1 --out-dir ${prefix}.spades_out --threads ${task.cpus} --scaffold
    else
      #canu -d ${prefix}.spades_out -p ${prefix} genomeSize=4m corOutCoverage=999 corMhapSensitivity=high useGrid=false maxThreads=${task.cpus} maxMemory=${task.memory.toGiga()}g -nanopore $R1
      flye --nano-raw $R1 --out-dir ${prefix}.spades_out --threads ${task.cpus} --scaffold --meta
    fi
    #ln -s ${prefix}.spades_out/${prefix}.contigs.fasta ${prefix}.contigs.fasta # for canu
    cut -f1,2,3 ${prefix}.spades_out/assembly_info.txt | awk -F'\t' 'NR>1{print \$1"\t"\$1"_length_"\$2"_cov_"\$3}' > flyeID_spadesID.txt
    fasta_tool --swap_ids flyeID_spadesID.txt ${prefix}.spades_out/assembly.fasta > ${prefix}.contigs.fasta
    ##ln -s ${prefix}.spades_out/assembly.fasta ${prefix}.contigs.fasta # for flye
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/ len=.*\$//g' > ${prefix}.ctgs.fasta
    """
}

} else {

process spades {
    tag "${prefix}"
    publishDir path: "${params.outdir}/spades", mode: 'copy'

    input:
    file clean_reads from normalized_reads_for_assembly

    output:
    file "${prefix}.ctg200.fasta" into contigs_for_nt, contigs_for_split
    file "${prefix}.ctgs.fasta" into contigs_for_quast1, contigs_for_quast2, contigs_for_checkm, contigs_for_prokka, contigs_for_prodigal, contigs_for_resfinder, contigs_for_pointfinder, contigs_for_tsne, contigs_for_augustus, contigs_for_eukcc

    when:
    params.ass

    script:
    prefix = clean_reads[0].toString() - ~/(_trimmed)?(_norm)?(_combined)?(\.R1)?(_1)?(_R1)?(\.1_val_1)?(_1_val_1)?(_val_1)?(_R1_val_1)?(\.fq)?(\.fastq)?(\.gz)?(\.bz2)?$/
    R1 = clean_reads[0].toString()
    mode = params.bulk ? "bulk" : "mda"
    if (single_end) {
    """
    if [ \"${mode}\" == \"bulk\" ]; then
      spades.py -s $R1 --careful --cov-cutoff auto -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    else
      spades.py --sc -s $R1 --careful -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    fi
    ln -s ${prefix}.spades_out/contigs.fasta ${prefix}.contigs.fasta
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/_length.*\$//g' > ${prefix}.ctgs.fasta
    """
    } else {
    R2 = clean_reads[1].toString()
    """
    if [ \"${mode}\" == \"bulk\" ]; then
      spades.py -1 $R1 -2 $R2 --careful --cov-cutoff auto -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    else
      spades.py --sc -1 $R1 -2 $R2 --careful -t ${task.cpus} -m ${task.memory.toGiga()} -o ${prefix}.spades_out
    fi
    ln -s ${prefix}.spades_out/contigs.fasta ${prefix}.contigs.fasta
    faFilterByLen.pl ${prefix}.contigs.fasta 200 > ${prefix}.ctg200.fasta
    cat ${prefix}.ctg200.fasta | sed 's/_length.*\$//g' > ${prefix}.ctgs.fasta
    """
    }
}

}

/*
 * STEP 7 - Evaluation using QUAST
 */
process quast_ref {
    label "quast"
    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    file fasta from fasta
    file gff from gff
    file ("*") from contigs_for_quast1.collect()
    file ("*") from bam_for_quast.collect()
    file ("*") from bai_for_quast.collect()

    output:
    file "quast/report.tsv" into quast_report1
    file "quast"

    when:
    denovo == false

    script:
    euk_cmd = euk ? ( params.fungus ? "--fungus" : "-e") : ""
    ref = fasta.exists() ? "-r $fasta" : ""
    gene = gff.exists() ? "--features gene:$gff" : ""
    """
    contigs=\$(ls *.ctgs.fasta | paste -sd " " -)
    labels=\$(ls *.ctgs.fasta | paste -sd "," - | sed 's/.ctgs.fasta//g')
    bams=\$(ls *.markdup.bam | paste -sd "," -)
    quast.py -o quast $ref $gene -m 200 -t ${task.cpus} $euk_cmd --rna-finding --bam \$bams -l \$labels --no-sv --no-read-stats \$contigs
    """
}

process quast_denovo {
    label "quast"
    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    file ("*") from contigs_for_quast2.collect()

    output:
    file "quast/report.tsv" into quast_report2
    file "quast"

    when:
    denovo == true

    script:
    euk_cmd = euk ? ( params.fungus ? "--fungus" : "-e") : ""
    """
    contigs=\$(ls *.ctgs.fasta | paste -sd " " -)
    labels=\$(ls *.ctgs.fasta | paste -sd "," - | sed 's/.ctgs.fasta//g')
    quast.py -o quast -m 200 -t ${task.cpus} $euk_cmd --rna-finding -l \$labels --no-sv --no-read-stats \$contigs
    """
}

/*
 * STEP 8.1 - Completeness and contamination evaluation using CheckM
 */
process checkm {
   publishDir "${params.outdir}/CheckM", mode: 'copy'

   input:
   file ('spades/*') from contigs_for_checkm.collect()

   output:
   file 'spades_checkM.txt'

   when:
   !euk

   script:
   checkm_wf = params.genus ? "taxonomy_wf" : "lineage_wf"
   """
   if [ \"${checkm_wf}\" == \"taxonomy_wf\" ]; then
     checkm taxonomy_wf -t ${task.cpus} --tab_table -f spades_checkM.txt -x fasta genus ${params.genus} spades spades_checkM
   else
     checkm lineage_wf -t ${task.cpus} -r --tab_table -f spades_checkM.txt -x fasta spades spades_checkM
   fi
   """
}

/*
 * STEP 9.1 - Annotate contigs using NT database
 */
process blast_nt {
   tag "${prefix}"
   publishDir "${params.outdir}/blob", mode: 'copy'

   input:
   file contigs from contigs_for_nt
   file nt from nt_db

   output:
   file "${prefix}_nt.out" into nt_for_blobtools_original
   file "${contigs}" into contigs_for_uniprot

   when:
   params.nt_db

   script:
   prefix = contigs.toString() - ~/(\.ctg200\.fasta)?(\.ctg200)?(\.fasta)?(\.fa)?$/
   """
   export BLASTDB=$nt
   blastn -query $contigs -db $nt/nt -outfmt '6 qseqid staxids bitscore std' \
     -max_target_seqs 1 -max_hsps 1 -evalue ${params.evalue} \
     -num_threads ${task.cpus} -out ${prefix}_nt.out
   """
}

/*
 * STEP 9.2 - Annotate contigs using Uniprot proteomes database
 */
process diamond_uniprot {
   tag "${prefix}"
   publishDir "${params.outdir}/blob", mode: 'copy'

   input:
   file contigs from contigs_for_uniprot
   file nt_out from nt_for_blobtools_original
   file uniprot from uniprot_db
   file "uniprot.taxids" from uniprot_taxids

   output:
   file "${prefix}_uniprot.taxified.out" into uniprot_for_blobtools
   file "${contigs}" into contigs_for_blob
   file "${nt_out}" into nt_for_blobtools
   val used into uniprot_real
   file "${prefix}_uniprot.*"

   script:
   prefix = contigs.toString() - ~/(\.ctg200\.fasta)?(\.ctg200)?(\.fasta)?(\.fa)?$/
   if ( uniprot.toString().equals("/dev/null") || uniprot.toString().equals("null") ) {
     used = false
     """
     touch ${prefix}_uniprot.out
     touch ${prefix}_uniprot.taxified.out
     """
   } else {
     used = true
     """
     diamond blastx --query $contigs --db $uniprot -p ${task.cpus} -o ${prefix}_uniprot.out \
       --outfmt 6 --sensitive --max-target-seqs 1 --evalue ${params.evalue} -b ${params.blockSize}
     blobtools taxify -f ${prefix}_uniprot.out -m uniprot.taxids -s 0 -t 2
     """
   }
}

/*
 * STEP 10.1 - Blobplot
 */
process blobtools {
   tag "${prefix}"
   publishDir "${params.outdir}/blob", mode: 'copy'

   input:
   file contigs from contigs_for_blob
   file anno from nt_for_blobtools
   val has_uniprot from uniprot_real
   file uniprot_anno from uniprot_for_blobtools

   output:
   file "${prefix}/${prefix}.blobDB*table.txt" into blob_tax
   file "${contigs}" into contigs_for_acdc
   file "${prefix}" into blob_tax_for_split

   script:
   prefix = contigs.toString() - ~/(\.ctg200\.fasta)?(\.ctg200)?(\.fasta)?(\.fa)?$/
   uniprot_anno_cmd = has_uniprot ? "-t $uniprot_anno" : ""
   """
   mkdir -p ${prefix}
   blobtools create -i $contigs -y spades -t $anno $uniprot_anno_cmd -o ${prefix}/${prefix} \
     --db /opt/conda/envs/nf-core-gongyh-scgs/lib/python3.6/site-packages/data/nodesDB.txt
   blobtools view -i ${prefix}/${prefix}.blobDB.json -r all -o ${prefix}/
   blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r phylum --format pdf -o ${prefix}/
   blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r order --format pdf -o ${prefix}/
   blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r family --format pdf -o ${prefix}/
   blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r genus --format pdf -o ${prefix}/
   blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r species --format pdf -o ${prefix}/
   blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r phylum --format png -o ${prefix}/
   blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r order --format png -o ${prefix}/
   blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r family --format png -o ${prefix}/
   blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r genus --format png -o ${prefix}/
   blobtools plot -i ${prefix}/${prefix}.blobDB.json --filelabel --notitle -l 200 -r species --format png -o ${prefix}/
   """
}

/*
 * STEP 10.2 - ACDC
 */
process acdc {
    tag "${prefix}"
    publishDir "${params.outdir}/acdc", mode: 'copy'

    input:
    file contigs from contigs_for_acdc
    file db from kraken_db
    file tax from blob_tax

    output:
    file "${prefix}"

    when:
    false

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    """
    cat $tax | grep -v '^#' | cut -f1,18 > genus.txt
    /usr/local/bin/acdc -i $contigs -m 1000 -b 100 -o $prefix -K $db -x genus.txt -T ${task.cpus}
    """
}

/*
 * STEP 10.3 - tSNE
 */
process tsne {
    tag "${prefix}"
    publishDir "${params.outdir}/tsne", mode: 'copy'

    input:
    file contigs from contigs_for_tsne

    output:
    file "${prefix}_tsne.tsv"

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    """
    faFilterByLen.pl ${contigs} 1000 > ${prefix}.ctg1k.fasta
    if [ -s ${prefix}.ctg1k.fasta ]
    then
      kpal count -k 4 -r ${prefix}.ctg1k.fasta ${prefix}.4mer
      kmer_tsne.py ${prefix}.4mer ${prefix}_tsne.tsv ${task.cpus}
    else
      touch ${prefix}_tsne.tsv
    fi
    """
}

if (!euk) {
/*
 * STEP 11 - Find genes using Prokka
 */
process prokka {
   tag "$prefix"
   publishDir "${params.outdir}/prokka", mode: 'copy'

   input:
   file contigs from contigs_for_prokka

   output:
   file "$prefix" into prokka_for_mqc1, prokka_for_mqc2, prokka_for_split
   file "$prefix/${prefix}.faa" into faa_eggnog, faa_kofam

   when:
   !euk

   script:
   prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
   """
   if [ \$(id -u) -eq 0 ]; then
     wget -c -q -t 1 -T 60 ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz -O linux64.tbl2asn.gz && gunzip linux64.tbl2asn.gz && chmod +x linux64.tbl2asn && mv linux64.tbl2asn /opt/conda/envs/nf-core-gongyh-scgs/bin/tbl2asn
   fi
   cat $contigs | sed 's/_length.*\$//g' > ${prefix}_node.fa
   prokka --outdir $prefix --prefix $prefix --addgenes --cpus ${task.cpus} ${prefix}_node.fa || echo "Ignore minor errors of prokka!"
   sed '/^##FASTA/Q' ${prefix}/${prefix}.gff > ${prefix}/${prefix}_noseq.gff
   gff2bed < ${prefix}/${prefix}_noseq.gff | cut -f1,4 | grep -v gene > ${prefix}/${prefix}_ctg_genes.tsv
   prokka_postprocess.py ${prefix}/${prefix}_ctg_genes.tsv ${prefix}/${prefix}.tsv > ${prefix}/${prefix}_all.tsv
   """
}
/*
 * STEP 11.2 - Find genes using prodigal
 */
process prodigal {
   tag "$prefix"
   publishDir "${params.outdir}/prodigal", mode: 'copy'

   input:
   file contigs from contigs_for_prodigal

   output:
   file "$prefix"

   when:
   !euk

   script:
   prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
   """
   mkdir -p ${prefix}
   prodigal -i $contigs -o ${prefix}/${prefix}.gbk -a ${prefix}/${prefix}.proteins.faa -p meta
   """
}


} else {
prokka_for_mqc1 = file('/dev/null')
prokka_for_mqc2 = file('/dev/null')
prokka_for_split = file('/dev/null')
/*
 * STEP 11.2 - Find genes using Augustus
 */
process augustus {
   tag "$prefix"
   publishDir "${params.outdir}/augustus", mode: 'copy'

   input:
   file contigs from contigs_for_augustus

   output:
   file "${prefix}.aa" into faa_eukcc, faa_eggnog, faa_kofam
   file "${prefix}*"

   when:
   euk

   script:
   prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
   """
   # clean id
   cat $contigs | sed 's/_length.*\$//g' > ${prefix}_clean.fasta
   # mask genome
   tantan ${prefix}_clean.fasta > ${prefix}_mask.fasta
   # gene prediction
   augustus --species=${params.augustus_species} --gff3=on --uniqueGeneId=true --protein=on --codingseq=on ${prefix}_mask.fasta > ${prefix}.gff
   # generate proteins
   getAnnoFasta.pl ${prefix}.gff
   """
}

/*
 * STEP 11.3 - Completeness and contamination evaluation using EukCC for euk
 */
process eukcc {
   publishDir "${params.outdir}/EukCC", mode: 'copy'

   input:
   file faa from faa_eukcc
   file db from eukcc_db

   output:
   file "${prefix}"

   when:
   euk && eukcc_db

   script:
   prefix = faa.toString() - ~/(\.faa)?(\.aa)?(\.fasta)?(\.fa)?$/
   """
   export HOME=/tmp/
   if [ -f "/tmp/.etetoolkit/taxa.sqlite" ]; then
     echo "NCBI taxa database exist!"
   else
     python -c "from ete3 import NCBITaxa; ncbi = NCBITaxa(taxdump_file='/opt/nf-core-scgs/taxdump.tar.gz')"
   fi
   eukcc --db ${db} --ncores ${task.cpus} --outdir ${prefix} --protein ${faa} || echo "Ignore minor errors of eukcc!"
   """
}

}


/*
 * STEP 12 - MultiQC
 */
preseq_for_multiqc = file('/dev/null')
process multiqc_ref {
    label "multiqc"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_multiqc_config1
    file ('fastqc/*') from fastqc_results1.collect()
    file ('software_versions/*') from software_versions_yaml1
    file ('trimgalore/*') from trimgalore_results1.collect()
    file ('fastqc2/*') from trimgalore_fastqc_reports1.collect()
    file ('samtools/*') from samtools_stats.collect()
    file ('preseq/*') from preseq_for_multiqc.collect()
    file ('*') from qualimap_for_multiqc
    file ('quast/*') from quast_report1
    file ('prokka/*') from prokka_for_mqc1.collect()
    file ('kraken/*') from kraken_for_mqc1.collect()
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report1
    file "*_data"

    when:
    denovo == false

    script:
    """
    multiqc -f --config $multiqc_config .
    """
}

process multiqc_denovo {
    label "multiqc"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_multiqc_config2
    file ('fastqc/*') from fastqc_results2.collect()
    file ('software_versions/*') from software_versions_yaml2
    file ('trimgalore/*') from trimgalore_results2.collect()
    file ('fastqc2/*') from trimgalore_fastqc_reports2.collect()
    file ('quast/*') from quast_report2.ifEmpty('/dev/null')
    file ('prokka/*') from prokka_for_mqc2.collect()
    file ('kraken/*') from kraken_for_mqc2.collect()
    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into multiqc_report2
    file "*_data"

    when:
    denovo == true

    script:
    """
    multiqc -f --config $multiqc_config .
    """
}

/*
 * STEP 13 - Annotate genes using EggNOG
 */
process eggnog {
   tag "$prefix"
   publishDir "${params.outdir}/eggnog", mode: 'copy'

   input:
   file faa from faa_eggnog
   file db from eggnog_db

   output:
   file "${prefix}.emapper.annotations"

   when:
   eggnog_db

   script:
   prefix = faa.toString() - ~/(\.proteins\.fa)?(\.faa)?$/
   """
   set +u
   source activate base
   emapper.py -i $faa -o $prefix --data_dir $db --dmnd_db $db/eggnog_proteins.dmnd -m diamond --cpu ${task.cpus}
   """
}

/*
 * STEP 13.1 - Annotate genes using KOfamKOALA
 */
process kofam {
   tag "$prefix"
   publishDir "${params.outdir}/kofam", mode: 'copy'

   input:
   file faa from faa_kofam
   file profile from kofam_profile
   file ko_list from kofam_kolist

   output:
   file "${prefix}_KOs_*.txt" into kofam_for_split

   when:
   kofam_profile && kofam_kolist

   script:
   prefix = faa.toString() - ~/(\.proteins\.fa)?(\.faa)?$/
   """
   exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -o ${prefix}_KOs_detail.txt ${faa}
   exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -r -f mapper -o ${prefix}_KOs_mapper.txt ${faa}
   exec_annotation -p ${profile} -k ${ko_list} --cpu ${task.cpus} -T 0.8 --keep-tabular -r -f mapper-one-line -o ${prefix}_KOs_mapper2.txt ${faa}
   kofam_postprocess.py /opt/nf-core-scgs/assets/ko_KO.txt ${prefix}_KOs_mapper.txt > ${prefix}_KOs_ko.txt
   """
}

/*
 * STEP 14 - Find ARGs
 */
process resfinder {
    tag "$prefix"
    publishDir "${params.outdir}/ARG", mode: 'copy'

    input:
    file contigs from contigs_for_resfinder
    file db from resfinder_db

    output:
    file "${prefix}/*"

    when:
    !euk && params.acquired && resfinder_db

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    """
    mkdir -p $prefix
    python /opt/resfinder/resfinder.py -i $contigs -o $prefix -p $db -mp blastn -x
    rm -rf $prefix/tmp
    """
}

/*
 * STEP 15 - Find point mutations
 */
process pointfinder {
    tag "$prefix"
    publishDir "${params.outdir}/ARG", mode: 'copy'

    input:
    file contigs from contigs_for_pointfinder
    file db from pointfinder_db

    output:
    file "${prefix}/*"

    when:
    !euk && params.point && pointfinder_db

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    species = params.pointfinder_species
    known_snp = params.only_known ? "" : "-l 0.4 -r all -u"
    """
    mkdir -p $prefix
    python /opt/pointfinder/PointFinder.py -p $db \
      -m blastn -m_p /opt/conda/bin/blastn $known_snp \
      -i $contigs -o $prefix -s $species
    rm -rf $prefix/tmp
    """
}


/*
 * Split the genome by contig annotation and checkm for each genus
 */

process split_checkm {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    file "results/spades/*" from contigs_for_split.collect()
    file "results/blob/*" from blob_tax_for_split.collect()
    file "results/prokka/*" from prokka_for_split.collect()
    file "results/kofam/*" from kofam_for_split.collect()

    output:
    file "split/*"

    when:
    params.split

    script:
    """
    cli.py tools scgs_split
    cd split
    samples=(`ls -d *_genus | sed 's/_genus//g'`)
    for sample in \${samples[*]}; do
      mkdir -p \${sample}_genus_checkM
      checkm lineage_wf -t ${task.cpus} -f \${sample}_genus_checkM.txt -x fasta \${sample}_genus \${sample}_genus_checkM || echo "Ignore internal errors!" 
    done
    """
}

/*
 * STEP 16 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py -o results_description.html $output_docs
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[gongyh/nf-core-scgs] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[gongyh/nf-core-scgs] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.maxMultiqcEmailFileSize)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            if (multiqc_report1 != null) {
                mqc_report = multiqc_report1.getVal()
            } else if (multiqc_report2 != null) {
                mqc_report = multiqc_report2.getVal()
            }
            if (mqc_report.getClass() == ArrayList){
                log.warn "[gongyh/nf-core-scgs] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[gongyh/nf-core-scgs] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[gongyh/nf-core-scgs] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[gongyh/nf-core-scgs] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    if(workflow.success){
        log.info "${c_purple}[gongyh/nf-core-scgs]${c_green} Pipeline complete${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[gongyh/nf-core-scgs]${c_red} Pipeline completed with errors${c_reset}"
    }

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
