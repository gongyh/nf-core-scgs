include { METABULI_TAXA } from '../../modules/local/metabuli_taxa'
include { METABULI_TAXONOMY_TABLES } from '../../modules/local/metabuli_taxonomy_tables'

workflow METABULI_CONTIG_TAXONOMY {
    take:
    contigs
    database

    main:
    classifications = METABULI_TAXA(contigs, database)
    tables = METABULI_TAXONOMY_TABLES(classifications.map { result -> tuple(result.meta, result.classifications) })
    ch_published = classifications.map { result -> [destination: 'taxonomy/metabuli', files: result] }
        .mix(tables.map { result -> [destination: 'taxonomy/metabuli', files: result] })

    emit:
    taxonomy = tables.map { result -> tuple(result.meta, result.taxonomy) }
    semibin_taxonomy = tables.map { result -> result.semibin_taxonomy }
    mqc_tsv = tables.map { result -> result.mqc_tsv }
    published = ch_published
}
