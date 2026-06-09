include { METABULI_TAXA } from '../../modules/local/metabuli_taxa'
include { VAMB_BIN } from '../../modules/local/taxvamb'
include { CLUSTERS_TO_SCAFFOLDS2BIN } from '../../modules/local/clusters_to_scaffolds2bin'

workflow TAXVAMB_INTEGRATION {
    take:
    ch_assembly
    ch_bams_stream
    ch_abundance

    main:
    if (!params.metabuli_db) {
        emit:
        scaffolds2bin = Channel.empty()
        mqc_tsv = Channel.empty()
        versions = Channel.empty()
        return
    }

    ch_assembly_single = ch_assembly.collect()
    ch_abundance_single = ch_abundance.collect()
    ch_bams_list = ch_bams_stream.collect()
    def meta = [id: 'merged']
    ch_assembly_tuple = ch_assembly.map { asm -> [meta, asm] }
    if (!params.metabuli_db) {
        error "METABULI_TAXA requires a database path (--metabuli_db). Please provide it."
    }
    ch_taxonomy = METABULI_TAXA(ch_assembly_tuple, file(params.metabuli_db, type: 'dir')).taxonomy

    ch_taxonomy_path = ch_taxonomy.map { _meta, tax -> tax }
    ch_bams_safe = ch_bams_list.map { bams -> [ bams ] }
    ch_vamb_input = ch_assembly_single
        .combine(ch_abundance_single)
        .combine(ch_bams_safe)
        .combine(ch_taxonomy_path)
        .map { row ->
            def vamb_meta = [id: 'merged']
            return [ vamb_meta, row[0], row[1], row[2], row[3] ]
        }
    VAMB_BIN( ch_vamb_input )
    ch_scaffolds2bin = VAMB_BIN.out.scaffolds2bin
    ch_cluster_file = VAMB_BIN.out.clusters_unsplit.map { _meta, file -> file }
    CLUSTERS_TO_SCAFFOLDS2BIN( ch_cluster_file )

    emit:
    scaffolds2bin = VAMB_BIN.out.scaffolds2bin
    mqc_tsv       = Channel.empty()
    versions      = VAMB_BIN.out.versions_vamb
}
