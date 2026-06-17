def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the database preparation pipeline is as follows:

    nextflow run gongyh/nf-core-scgs --prepare_databases -profile docker

    Options:
    --outdir                      The output directory where the databases will be saved (Default: ./databases)
    --db_type                     Comma-separated list of databases to prepare. Options: mmseqs, checkm2, kofam, eggnog, kraken2, gtdb, blob, metabuli, genomad, nt, all

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

    // MMseqs2 database
    if (db_types.contains("all") || db_types.contains("mmseqs")) {
        MMSEQS_DOWNLOAD()
        ch_versions = ch_versions.mix(MMSEQS_DOWNLOAD.out.versions)
        log.info "Prepared MMseqs2 database: ${params.outdir}/mmseqs_db"
    }

    // CheckM2 database
    if (db_types.contains("all") || db_types.contains("checkm2")) {
        CHECKM2_DOWNLOAD()
        ch_versions = ch_versions.mix(CHECKM2_DOWNLOAD.out.versions)
        log.info "Prepared CheckM2 database: ${params.outdir}/checkm2_db"
    }

    // KOfam database
    if (db_types.contains("all") || db_types.contains("kofam")) {
        KOFAM_DOWNLOAD()
        ch_versions = ch_versions.mix(KOFAM_DOWNLOAD.out.versions)
        log.info "Prepared KOfam database: ${params.outdir}/kofam_db"
    }

    // EggNOG database
    if (db_types.contains("all") || db_types.contains("eggnog")) {
        EGGNOG_DOWNLOAD()
        ch_versions = ch_versions.mix(EGGNOG_DOWNLOAD.out.versions)
        log.info "Prepared EggNOG database: ${params.outdir}/eggnog_db"
    }

    // Kraken2 database
    if (db_types.contains("all") || db_types.contains("kraken2")) {
        KRAKEN2_DOWNLOAD()
        ch_versions = ch_versions.mix(KRAKEN2_DOWNLOAD.out.versions)
        log.info "Prepared Kraken2 database: ${params.outdir}/kraken2_db"
    }

    // GTDB database
    if (db_types.contains("all") || db_types.contains("gtdb")) {
        GTDB_DOWNLOAD()
        ch_versions = ch_versions.mix(GTDB_DOWNLOAD.out.versions)
        log.info "Prepared GTDB database: ${params.outdir}/gtdb_db"
    }

    // Blobtools database
    if (db_types.contains("all") || db_types.contains("blob")) {
        BLOB_DOWNLOAD()
        ch_versions = ch_versions.mix(BLOB_DOWNLOAD.out.versions)
        log.info "Prepared Blobtools database: ${params.outdir}/blob_db"
    }

    // MetaBuli database
    if (db_types.contains("all") || db_types.contains("metabuli")) {
        METABULI_DOWNLOAD()
        ch_versions = ch_versions.mix(METABULI_DOWNLOAD.out.versions)
        log.info "Prepared MetaBuli database: ${params.outdir}/metabuli_db"
    }

    // GENOMAD database
    if (db_types.contains("all") || db_types.contains("genomad")) {
        GENOMAD_DOWNLOAD()
        ch_versions = ch_versions.mix(GENOMAD_DOWNLOAD.out.versions)
        log.info "Prepared GENOMAD database: ${params.outdir}/genomad_db"
    }

    // NCBI nt database
    if (db_types.contains("all") || db_types.contains("nt")) {
        NT_DOWNLOAD()
        ch_versions = ch_versions.mix(NT_DOWNLOAD.out.versions)
        log.info "Prepared NCBI nt database: ${params.outdir}/nt_db"
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
