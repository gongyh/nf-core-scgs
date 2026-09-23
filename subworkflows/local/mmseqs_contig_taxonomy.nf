include { MMSEQS_CREATEDB  } from '../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_TAXONOMY  } from '../../modules/nf-core/mmseqs/taxonomy/main'
include { MMSEQS_CREATETSV } from '../../modules/nf-core/mmseqs/createtsv/main'
include { MMSEQS_TAXONOMY_MULTIQC } from '../../modules/local/mmseqs_taxonomy_multiqc'

workflow MMSEQS_CONTIG_TAXONOMY {

    take:
    contigs            // channel: tuple val(meta), path(contigs)
    mmseqs_databases   // channel: path(mmseqs2 local db)

    main:
    ch_mmseqs_db              = channel.empty()
    ch_taxonomy_querydb       = channel.empty()
    ch_taxonomy_querydb_taxdb = channel.empty()
    ch_taxonomy_tsv           = channel.empty()
    ch_published              = channel.empty()

    // MMSEQS_DATABASE
    if ( mmseqs_databases != null ) {
        ch_mmseqs_db = mmseqs_databases
    } else {
        ch_mmseqs_db = channel.empty()
    }

    // Create db for query contigs, assign taxonomy and convert to table format
    // MMSEQS_CREATEDB
    MMSEQS_CREATEDB ( contigs )
    ch_taxonomy_querydb = MMSEQS_CREATEDB.out.db
    ch_published = ch_published.mix(ch_taxonomy_querydb.map { result -> [destination: 'binning/mmseqs2_taxa', files: result] })

    // MMSEQS_TAXONOMY
    MMSEQS_TAXONOMY ( ch_taxonomy_querydb, ch_mmseqs_db )
    ch_taxonomy_querydb_taxdb = MMSEQS_TAXONOMY.out.db_taxonomy
    ch_published = ch_published.mix(ch_taxonomy_querydb_taxdb.map { result -> [destination: 'binning/mmseqs2_taxa', files: result] })

    // MMSEQS_CREATETSV
    MMSEQS_CREATETSV ( ch_taxonomy_querydb_taxdb, [[:],[]], ch_taxonomy_querydb )
    ch_taxonomy_tsv = MMSEQS_CREATETSV.out.tsv
    MMSEQS_TAXONOMY_MULTIQC ( ch_taxonomy_tsv )
    ch_published = ch_published.mix(ch_taxonomy_tsv.map { result -> [destination: 'binning/mmseqs2_taxa', files: result] })

    emit:
    taxonomy    = ch_taxonomy_tsv           // channel: [ val(meta), tsv ]
    mqc_tsv     = MMSEQS_TAXONOMY_MULTIQC.out.map { result -> tuple(result.meta, result.mqc_tsv) }
    db_mmseqs   = ch_mmseqs_db              // channel: [ val(meta), mmseqs_database ]
    db_taxonomy = ch_taxonomy_querydb_taxdb // channel: [ val(meta), db_taxonomy ]
    db_contig   = ch_taxonomy_querydb       // channel: [ val(meta), db ]
    published   = ch_published
}
