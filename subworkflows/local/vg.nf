/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { VG_CONSTRUCT          } from '../../modules/local/vg/vg_construct'
include { VG_INDEX              } from '../../modules/local/vg/vg_index'
include { VG_CALL               } from '../../modules/local/vg/vg_call'

workflow VG {
    take:
    fasta
    trimmed_reads
    vcf

    main:
    ch_versions = Channel.empty()
    VG_CONSTRUCT (
        fasta,
        vcf
    )
    VG_INDEX (
        VG_CONSTRUCT.out.vg,
        trimmed_reads
    )
    VG_CALL (
        VG_CONSTRUCT.out.vg,
        VG_INDEX.out.gam
    )
    ch_versions = ch_versions.mix(VG_CALL.out.versions)

    emit:
    ch_versions
}
