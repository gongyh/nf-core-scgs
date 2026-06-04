include { MAKE_DUMMY_TAXONOMY } from '../../modules/local/make_dummy_taxonomy'
include { VAMB_BIN } from '../../modules/local/taxvamb'
include { PROCESS_VAMB_BINS } from '../../modules/local/process_vamb_bins'

workflow TAXVAMB_INTEGRATION {
    take:
    ch_assembly
    ch_bams_stream

    main:
    ch_bams_list = ch_bams_stream
        .toList()
        .map { bam_paths ->
            return bam_paths as List
        }
    //DUMMY
    ch_taxonomy = MAKE_DUMMY_TAXONOMY(ch_assembly).taxonomy

    // VAMB_BIN
    def vamb_meta = [id: 'merged']
    ch_vamb_input = ch_assembly
        .combine(ch_bams_list)
        .combine(ch_taxonomy)
        .map { assembly, bams, taxonomy ->
            [vamb_meta, assembly, [], bams, taxonomy]
        }
    VAMB_BIN( ch_vamb_input )
    ch_cluster_file = VAMB_BIN.out.clusters_unsplit.map { meta, file -> file }
    PROCESS_VAMB_BINS( ch_cluster_file )

    emit:
    scaffolds2bin = PROCESS_VAMB_BINS.out.scaffolds2bin
    mqc_tsv       = PROCESS_VAMB_BINS.out.mqc_tsv
    versions      = VAMB_BIN.out.versions_vamb.mix(PROCESS_VAMB_BINS.out.versions)
}
