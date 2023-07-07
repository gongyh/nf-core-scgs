/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { MINIMAP2_ALIGN        } from '../../modules/nf-core/minimap2/align/main'
include { BOWTIE2_ALIGN         } from '../../modules/nf-core/bowtie2/align/main'

workflow ALIGN {
    take:
    is_nanopore
    trimmed_reads
    bowtie2_index
    fasta

    main:
    ch_versions = Channel.empty()
    if ( is_nanopore ) {
        MINIMAP2_ALIGN (
            trimmed_reads,
            fasta,
            true,
            false,
            false
        )
        MINIMAP2_ALIGN.out.bam.set{bb_bam}
        ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    } else {
        BOWTIE2_ALIGN (
            trimmed_reads,
            bowtie2_index,
            false,
            true
        )
        BOWTIE2_ALIGN.out.bam.set{bb_bam}
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
    }

    emit:
    bb_bam
    ch_versions
}
