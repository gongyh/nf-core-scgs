include { METABULI_TAXA } from '../../modules/local/metabuli_taxa'
include { VAMB_BIN } from '../../modules/local/taxvamb'

workflow TAXVAMB_INTEGRATION {
    take:
    ch_assembly
    ch_coverage

    main:
    if (!params.metabuli_db) {
        emit:
        scaffolds2bin = Channel.empty()
        mqc_tsv = Channel.empty()
        versions = Channel.empty()
        return
    }

    ch_assembly_single = ch_assembly.collect()
    ch_coverage_single = ch_coverage.collect()
    def meta = [id: 'merged']
    ch_assembly_tuple = ch_assembly.map { asm -> [meta, asm] }
    if (!params.metabuli_db) {
        error "METABULI_TAXA requires a database path (--metabuli_db). Please provide it."
    }
    ch_taxonomy = METABULI_TAXA(ch_assembly_tuple, file(params.metabuli_db, type: 'dir')).taxonomy

    ch_taxonomy_path = ch_taxonomy.map { _meta, tax -> tax }
    ch_vamb_input = ch_assembly_single
        .combine(ch_coverage_single)
        .combine(ch_taxonomy_path)
        .map { row ->
            def vamb_meta = [id: 'merged']
            return [ vamb_meta, row[0], row[1], [], row[2] ]
        }
    VAMB_BIN( ch_vamb_input )

    emit:
    scaffolds2bin = VAMB_BIN.out.scaffolds2bin
    mqc_tsv       = Channel.empty()
    versions      = VAMB_BIN.out.versions_vamb
}
