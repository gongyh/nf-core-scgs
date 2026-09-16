nextflow.enable.types = true

include { METABULI_TAXA } from '../../modules/local/metabuli_taxa'
include { VAMB_BIN } from '../../modules/local/taxvamb'

workflow TAXVAMB_INTEGRATION {
    take:
    ch_assembly: Channel<Path>
    ch_coverage: Channel<Path>

    main:
    if (params.metabuli_db) {
        def meta = [id: 'merged']
        ch_assembly_tuple = ch_assembly.map { asm -> tuple(meta, asm) }
        metabuli_taxa = METABULI_TAXA(ch_assembly_tuple, file(params.metabuli_db, type: 'dir'))
        ch_taxonomy_path = metabuli_taxa.map { result -> result.taxonomy }
        ch_vamb_input = ch_assembly
            .combine(ch_coverage)
            .combine(ch_taxonomy_path)
            .map { assembly, coverage, taxonomy ->
                tuple([id: 'merged'], assembly, coverage, taxonomy)
            }
        vamb_bin = VAMB_BIN(ch_vamb_input)
        ch_scaffolds2bin = vamb_bin.map { result -> tuple(result.meta, result.scaffolds2bin) }
        ch_versions = vamb_bin.map { result -> result.versions }
    } else {
        ch_scaffolds2bin = channel.empty()
        ch_versions = channel.empty()
    }

    emit:
    scaffolds2bin: Channel<Tuple<Map,Path>> = ch_scaffolds2bin
    mqc_tsv: Channel<Path> = channel.empty()
    versions: Channel<Path> = ch_versions
}
