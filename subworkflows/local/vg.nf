nextflow.enable.types = true

/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { VG_CONSTRUCT          } from '../../modules/local/vg/vg_construct'
include { VG_INDEX              } from '../../modules/local/vg/vg_index'
include { VG_CALL               } from '../../modules/local/vg/vg_call'

workflow VG {
    take:
    fasta: Path
    trimmed_reads: Channel<Tuple<Map,List<Path>>>
    vcf: Path

    main:
    ch_versions = channel.empty()
    vg_construct = VG_CONSTRUCT(fasta, vcf)
    ch_index_input = trimmed_reads.combine(vg_construct.map { result -> result.vg })
    vg_index = VG_INDEX(ch_index_input)
    ch_call_input = vg_index
        .map { result -> tuple(result.meta, result.gam) }
        .combine(vg_construct.map { result -> result.vg })
    vg_call = VG_CALL(ch_call_input)
    ch_versions = ch_versions.mix(vg_call.map { result -> result.versions })

    emit:
    ch_versions: Channel<Path>
}
