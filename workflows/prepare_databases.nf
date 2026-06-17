def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the database preparation pipeline is as follows:

    nextflow run gongyh/nf-core-scgs --prepare_databases -profile docker

    Options:
    --outdir                      The output directory where the databases will be saved (Default: ./databases)
    --db_type                     Comma-separated list of databases to prepare. Options: mmseqs, checkm2, kofam, eggnog, kraken2, gtdb, blob, metabuli, genomad, nt, all

    Database download URLs:
    --mmseqs_db_url               URL for MMseqs2 database (Default: https://mmseqs.com/databases)
    --checkm2_db_url              URL for CheckM2 database (Default: https://data.ace.uq.edu.au/public/CheckM2/)
    --kofam_db_url                URL for KOfam database (Default: https://www.genome.jp/ftp/db/kofam/)
    --eggnog_db_url               URL for EggNOG database (Default: https://eggnog5.embl.de/download/eggnog_5.0/)
    --kraken2_db_url              URL for Kraken2 database (Default: https://genome-idx.s3.amazonaws.com/kraken)
    --gtdb_db_url                 URL for GTDB database (Default: https://data.gtdb.ecogenomic.org/releases/)
    --blob_db_url                 URL for Blobtools nodesDB (Default: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)
    --metabuli_db_url             URL for MetaBuli database (Default: https://github.com/khyox/metabuli)
    --genomad_db_url              URL for GENOMAD database (Default: https://portal.nersc.gov/GENOMAD/)
    --nt_db_url                   URL for NCBI nt database (Default: https://ftp.ncbi.nlm.nih.gov/blast/db/)

    Generic options:
    --help                        Display this help message
    --monochrome_logs             Do not use coloured log outputs
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

if (params.help){
    helpMessage()
    exit 0
}

// Initialize parameters with default values
params.outdir = "./databases"
params.db_type = "all"
params.mmseqs_db_url = "https://mmseqs.com/databases"
params.checkm2_db_url = "https://data.ace.uq.edu.au/public/CheckM2/"
params.kofam_db_url = "https://www.genome.jp/ftp/db/kofam/"
params.eggnog_db_url = "https://eggnog5.embl.de/download/eggnog_5.0/"
params.kraken2_db_url = "https://genome-idx.s3.amazonaws.com/kraken"
params.gtdb_db_url = "https://data.gtdb.ecogenomic.org/releases/"
params.blob_db_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
params.metabuli_db_url = "https://github.com/khyox/metabuli"
params.genomad_db_url = "https://portal.nersc.gov/GENOMAD/"
params.nt_db_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/"

/*
 * Import modules
 */
include { MMSEQS_DOWNLOAD         } from '../modules/local/mmseqs_download'
include { CHECKM2_DOWNLOAD        } from '../modules/local/checkm2_download'
include { KOFAM_DOWNLOAD          } from '../modules/local/kofam_download'
include { EGGNOG_DOWNLOAD         } from '../modules/local/eggnog_download'
include { KRAKEN2_DOWNLOAD        } from '../modules/local/kraken2_download'
include { GTDB_DOWNLOAD           } from '../modules/local/gtdb_download'
include { BLOB_DOWNLOAD           } from '../modules/local/blob_download'
include { METABULI_DOWNLOAD       } from '../modules/local/metabuli_download'
include { GENOMAD_DOWNLOAD        } from '../modules/local/genomad_download'
include { NT_DOWNLOAD             } from '../modules/local/nt_download'
include { GET_SOFTWARE_VERSIONS   } from '../modules/local/get_software_versions/main'

/*
 * Workflow
 */
workflow PREPARE_DATABASES {
    ch_versions = Channel.empty()

    def db_types = params.db_type.toLowerCase().split(',').collect { it.trim() }

    // Create output directories
    def mmseqs_out = file("${params.outdir}/mmseqs")
    def checkm2_out = file("${params.outdir}/checkm2")
    def kofam_out = file("${params.outdir}/kofam")
    def eggnog_out = file("${params.outdir}/eggnog")
    def kraken2_out = file("${params.outdir}/kraken2")
    def gtdb_out = file("${params.outdir}/gtdb")
    def blob_out = file("${params.outdir}/blob")
    def metabuli_out = file("${params.outdir}/metabuli")
    def genomad_out = file("${params.outdir}/genomad")
    def nt_out = file("${params.outdir}/nt")

    // MMseqs2 database
    if (db_types.contains("all") || db_types.contains("mmseqs")) {
        MMSEQS_DOWNLOAD(params.mmseqs_db_url, mmseqs_out)
        ch_versions = ch_versions.mix(MMSEQS_DOWNLOAD.out.versions)
        log.info "Prepared MMseqs2 database: ${mmseqs_out}"
    }

    // CheckM2 database
    if (db_types.contains("all") || db_types.contains("checkm2")) {
        CHECKM2_DOWNLOAD(params.checkm2_db_url, checkm2_out)
        ch_versions = ch_versions.mix(CHECKM2_DOWNLOAD.out.versions)
        log.info "Prepared CheckM2 database: ${checkm2_out}"
    }

    // KOfam database
    if (db_types.contains("all") || db_types.contains("kofam")) {
        KOFAM_DOWNLOAD(params.kofam_db_url, kofam_out)
        ch_versions = ch_versions.mix(KOFAM_DOWNLOAD.out.versions)
        log.info "Prepared KOfam database: ${kofam_out}"
    }

    // EggNOG database
    if (db_types.contains("all") || db_types.contains("eggnog")) {
        EGGNOG_DOWNLOAD(params.eggnog_db_url, eggnog_out)
        ch_versions = ch_versions.mix(EGGNOG_DOWNLOAD.out.versions)
        log.info "Prepared EggNOG database: ${eggnog_out}"
    }

    // Kraken2 database
    if (db_types.contains("all") || db_types.contains("kraken2")) {
        KRAKEN2_DOWNLOAD(params.kraken2_db_url, kraken2_out)
        ch_versions = ch_versions.mix(KRAKEN2_DOWNLOAD.out.versions)
        log.info "Prepared Kraken2 database: ${kraken2_out}"
    }

    // GTDB database
    if (db_types.contains("all") || db_types.contains("gtdb")) {
        GTDB_DOWNLOAD(params.gtdb_db_url, gtdb_out)
        ch_versions = ch_versions.mix(GTDB_DOWNLOAD.out.versions)
        log.info "Prepared GTDB database: ${gtdb_out}"
    }

    // Blobtools database
    if (db_types.contains("all") || db_types.contains("blob")) {
        BLOB_DOWNLOAD(params.blob_db_url, blob_out)
        ch_versions = ch_versions.mix(BLOB_DOWNLOAD.out.versions)
        log.info "Prepared Blobtools database: ${blob_out}"
    }

    // MetaBuli database
    if (db_types.contains("all") || db_types.contains("metabuli")) {
        METABULI_DOWNLOAD(params.metabuli_db_url, metabuli_out)
        ch_versions = ch_versions.mix(METABULI_DOWNLOAD.out.versions)
        log.info "Prepared MetaBuli database: ${metabuli_out}"
    }

    // GENOMAD database
    if (db_types.contains("all") || db_types.contains("genomad")) {
        GENOMAD_DOWNLOAD(params.genomad_db_url, genomad_out)
        ch_versions = ch_versions.mix(GENOMAD_DOWNLOAD.out.versions)
        log.info "Prepared GENOMAD database: ${genomad_out}"
    }

    // NCBI nt database
    if (db_types.contains("all") || db_types.contains("nt")) {
        NT_DOWNLOAD(params.nt_db_url, nt_out)
        ch_versions = ch_versions.mix(NT_DOWNLOAD.out.versions)
        log.info "Prepared NCBI nt database: ${nt_out}"
    }

    // GET_SOFTWARE_VERSIONS
    if (!ch_versions.isEmpty()) {
        GET_SOFTWARE_VERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
    }

    log.info "Database preparation completed. All databases saved to: ${params.outdir}"
}

def nfcoreHeader(){
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
    ${c_purple}  gongyh/nf-core-scgs Database Preparation v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}
