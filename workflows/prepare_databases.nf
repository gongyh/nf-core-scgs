nextflow.enable.types = true

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    nextflow run gongyh/nf-core-scgs --prepare_databases -profile docker

    Workflow selection:
    --prepare_databases           Run the database preparation workflow

    Database options:
    --db_type <list>              Comma-separated databases to prepare (default: all)
                                    Available: mmseqs, checkm2, kofam, eggnog, kraken2,
                                    gtdb, blob, metabuli, genomad, nt, all

    Output and execution:
    --outdir <path>               Output directory for prepared databases (default: ./databases)
    --monochrome_logs             Disable coloured log output
    --help                        Display this help message
    -profile                      Configuration profile(s), for example: docker, singularity, conda
    """.stripIndent()
}

/*
 * Import modules
 */
include { MMSEQS_DBDOWNLOAD         } from '../modules/local/mmseqs_download'
include { CHECKM2_DBDOWNLOAD        } from '../modules/local/checkm2_download'
include { KOFAM_DBDOWNLOAD          } from '../modules/local/kofam_download'
include { EGGNOG_DBDOWNLOAD         } from '../modules/local/eggnog_download'
include { KRAKEN2_DBDOWNLOAD        } from '../modules/local/kraken2_download'
include { GTDB_DBDOWNLOAD           } from '../modules/local/gtdb_download'
include { BLOB_DBDOWNLOAD           } from '../modules/local/blob_download'
include { METABULI_DBDOWNLOAD       } from '../modules/local/metabuli_download'
include { GENOMAD_DBDOWNLOAD        } from '../modules/local/genomad_download'
include { NT_DBDOWNLOAD             } from '../modules/local/nt_download'
include { GET_SOFTWARE_VERSIONS     } from '../modules/local/get_software_versions/main'

/*
 * Workflow
 */
workflow PREPARE_DATABASES {
    main:
    params.outdir = "./databases"
    params.db_type = "all"
    ch_published = channel.empty()

    def db_types = params.db_type.toLowerCase().split(',').collect { db_type -> db_type.trim() }

    // MMseqs2 database
    if (db_types.contains("all") || db_types.contains("mmseqs")) {
        mmseqs_download = MMSEQS_DBDOWNLOAD()
        ch_published = ch_published.mix(mmseqs_download.map { result -> [destination: '.', files: result.db] })
        log.info "Prepared MMseqs2 database: ${params.outdir}/mmseqs_db"
    }

    // CheckM2 database
    if (db_types.contains("all") || db_types.contains("checkm2")) {
        checkm2_download = CHECKM2_DBDOWNLOAD()
        ch_published = ch_published.mix(checkm2_download.map { result -> [destination: '.', files: result.db] })
        log.info "Prepared CheckM2 database: ${params.outdir}/checkm2_db"
    }

    // KOfam database
    if (db_types.contains("all") || db_types.contains("kofam")) {
        kofam_download = KOFAM_DBDOWNLOAD()
        ch_published = ch_published.mix(kofam_download.map { result -> [destination: '.', files: result.db] })
        log.info "Prepared KOfam database: ${params.outdir}/kofam_db"
    }

    // EggNOG database
    if (db_types.contains("all") || db_types.contains("eggnog")) {
        eggnog_download = EGGNOG_DBDOWNLOAD()
        ch_published = ch_published.mix(eggnog_download.map { result -> [destination: '.', files: result.db] })
        log.info "Prepared EggNOG database: ${params.outdir}/eggnog_db"
    }

    // Kraken2 database
    if (db_types.contains("all") || db_types.contains("kraken2")) {
        kraken2_download = KRAKEN2_DBDOWNLOAD()
        ch_published = ch_published.mix(kraken2_download.map { result -> [destination: '.', files: result.db] })
        log.info "Prepared Kraken2 database: ${params.outdir}/kraken2_db"
    }

    // GTDB database
    if (db_types.contains("all") || db_types.contains("gtdb")) {
        gtdb_download = GTDB_DBDOWNLOAD()
        ch_published = ch_published.mix(gtdb_download.map { result -> [destination: '.', files: result.db] })
        log.info "Prepared GTDB database: ${params.outdir}/gtdb_db"
    }

    // Blobtools database
    if (db_types.contains("all") || db_types.contains("blob")) {
        blob_download = BLOB_DBDOWNLOAD()
        ch_published = ch_published.mix(blob_download.map { result -> [destination: '.', files: result.db] })
        log.info "Prepared Blobtools database: ${params.outdir}/blob_db"
    }

    // MetaBuli database
    if (db_types.contains("all") || db_types.contains("metabuli")) {
        metabuli_download = METABULI_DBDOWNLOAD()
        ch_published = ch_published.mix(metabuli_download.map { result -> [destination: '.', files: result.db] })
        log.info "Prepared MetaBuli database: ${params.outdir}/metabuli_db"
    }

    // GENOMAD database
    if (db_types.contains("all") || db_types.contains("genomad")) {
        genomad_download = GENOMAD_DBDOWNLOAD()
        ch_published = ch_published.mix(genomad_download.map { result -> [destination: '.', files: result.db] })
        log.info "Prepared GENOMAD database: ${params.outdir}/genomad_db"
    }

    // NCBI nt database
    if (db_types.contains("all") || db_types.contains("nt")) {
        nt_download = NT_DBDOWNLOAD()
        ch_published = ch_published.mix(nt_download.map { result -> [destination: '.', files: result.db] })
        log.info "Prepared NCBI nt database: ${params.outdir}/nt_db"
    }

    // GET_SOFTWARE_VERSIONS
    software_versions = GET_SOFTWARE_VERSIONS (
        channel.topic('local_versions')
            .unique()
            .collectFile(name: 'collated_versions.yml', newLine: true)
    )
    ch_published = ch_published.mix(software_versions.map { result -> [destination: 'pipeline_info', files: [result.yml, result.mqc_yml]] })

    log.info "Database preparation completed. All databases saved to: ${params.outdir}"

    emit:
    published = ch_published
}

def nfcoreHeader(){
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
    ${c_purple}  gongyh/nf-core-scgs Database Preparation v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}
