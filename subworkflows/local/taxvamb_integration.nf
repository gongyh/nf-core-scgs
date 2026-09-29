include { VAMB_BIN } from '../../modules/local/taxvamb'

workflow TAXVAMB_INTEGRATION {
    take:
    ch_assembly
    ch_coverage
    ch_taxonomy

    main:
    ch_published = channel.empty()
    ch_taxonomy_path = ch_taxonomy.map { _meta, taxonomy -> taxonomy }
    ch_vamb_input = ch_assembly
        .combine(ch_coverage)
        .combine(ch_taxonomy_path)
        .map { assembly, coverage, taxonomy ->
            tuple([id: 'merged'], assembly, coverage, taxonomy)
        }
    vamb_bin = VAMB_BIN(ch_vamb_input)
    ch_published = ch_published.mix(vamb_bin.map { result -> [destination: 'binning/taxvamb', files: result] })
    ch_scaffolds2bin = vamb_bin.map { result -> tuple(result.meta, result.scaffolds2bin) }
    ch_mqc_tsv = vamb_bin.map { result -> result.mqc_tsv }

    emit:
    scaffolds2bin = ch_scaffolds2bin
    mqc_tsv = ch_mqc_tsv
    published = ch_published
}
