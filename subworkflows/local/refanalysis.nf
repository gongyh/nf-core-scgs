/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { SAMTOOLS              } from '../../modules/local/samtools'
include { PRESEQ                } from '../../modules/local/preseq'
include { INDELREALIGN          } from '../../modules/local/indelrealign'
include { MONOVAR               } from '../../modules/local/monovar'
include { QUAST_REF             } from '../../modules/local/quast_ref'
include { ANEUFINDER            } from '../../modules/local/aneufinder'
include { CIRCLIZE              } from '../../modules/local/circlize'
include { SAVE_REFERENCE        } from '../../modules/local/save_reference'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { BOWTIE2_ALIGN         } from '../../modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN        } from '../../modules/nf-core/minimap2/align/main'
include { QUALIMAP_BAMQC        } from '../../modules/nf-core/qualimap/bamqc/main'

workflow REFANALYSIS {
    take:
    trimmed_reads
    fasta
    gff
    is_nanopore
    bowtie2_index
    is_vcf
    vcf
    is_snv
    is_cnv
    is_bulk
    is_single

    main:
    ch_versions = Channel.empty()
    ch_multiqc = Channel.empty()
    bb_bam = Channel.empty()
    SAVE_REFERENCE ( fasta, gff )
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

    if (is_vcf ) {
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
    }

    ch_multiqc_samtools = Channel.empty()
    ch_multiqc_preseq   = Channel.empty()
    ch_multiqc_qualimap = Channel.empty()
    quast_bam = Channel.empty()
    quast_bai = Channel.empty()
    SAMTOOLS (
        bb_bam,
        SAVE_REFERENCE.out.bed
    )
    SAMTOOLS.out.bam.set{quast_bam}
    SAMTOOLS.out.bai.set{quast_bai}
    ch_versions = ch_versions.mix(SAMTOOLS.out.versions)
    ch_multiqc_samtools = SAMTOOLS.out.stats
    ch_multiqc = ch_multiqc.mix(ch_multiqc_samtools.collect{it[1]}.ifEmpty([]))
    if (!is_nanopore) {
        PRESEQ ( SAMTOOLS.out.bed )
        ch_versions = ch_versions.mix(PRESEQ.out.versions)
        ch_multiqc_preseq = PRESEQ.out.txt
        ch_multiqc = ch_multiqc.mix(ch_multiqc_preseq.collect{it[1]}.ifEmpty([]))
    }

    QUALIMAP_BAMQC (
        SAMTOOLS.out.bam,
        gff
    )
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)
    ch_multiqc_qualimap = QUALIMAP_BAMQC.out.results
    ch_multiqc = ch_multiqc.mix(ch_multiqc_qualimap.collect{it[1]}.ifEmpty([]))
    if (is_snv && !is_nanopore) {
        INDELREALIGN (
            SAMTOOLS.out.bam,
            fasta
        )
        ch_versions = ch_versions.mix(INDELREALIGN.out.versions)
    }

    if (!is_bulk && is_snv && !is_nanopore) {
        MONOVAR (
            INDELREALIGN.out.bam.collect{it[1]},
            INDELREALIGN.out.bai.collect{it[1]},
            fasta
        )
        ch_versions = ch_versions.mix(MONOVAR.out.versions)
    }

    if (!is_bulk && is_cnv && !is_single && !is_nanopore) {
        ANEUFINDER (
            SAMTOOLS.out.bam.collect{it[1]},
            SAMTOOLS.out.bai.collect{it[1]}
        )
        ch_versions = ch_versions.mix(ANEUFINDER.out.versions)
    }

    CIRCLIZE (
        SAMTOOLS.out.bed,
        SAVE_REFERENCE.out.bed
    )
    ch_versions = ch_versions.mix(CIRCLIZE.out.versions)

    emit:
    ch_versions
    bam = quast_bam.collect{it[1]}
    bai = quast_bai.collect{it[1]}
    ch_multiqc
}
