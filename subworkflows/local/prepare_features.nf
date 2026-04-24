include { PANDEPTH_COVERAGE          } from '../../modules/local/pandepth_coverage'
include { PANDEPTH_MERGE             } from '../../modules/local/pandepth_merge'
include { PRODIGAL                   } from '../../modules/local/prodigal'
include { KMER_COUNT                 } from '../../modules/local/kmer_count'
include { SUMMARIZE_FEATURE_MATRIX   } from '../../modules/local/summarize_feature_matrix'

workflow PREPARE_FEATURES {
    take:
    ch_fasta
    ch_fai
    ch_bams

    main:
    ch_versions = Channel.empty()

    ch_fasta_file = ch_fasta
        .map { it -> it instanceof List ? it : [it] }
        .flatten()
        .filter { it.toString().endsWith('.fasta') || it.toString().endsWith('.fa') }
        .first()
    ch_fai_file = ch_fai
        .map { it -> it instanceof List ? it : [it] }
        .flatten()
        .filter { it.toString().endsWith('.fai') }
        .first()
    ch_bams_with_bai = ch_bams.map { meta, bam -> [meta, bam, []] }
    ch_depth = PANDEPTH_COVERAGE( ch_bams_with_bai, ch_fasta_file, ch_fai_file ).depth
    ch_all_depth = ch_depth.map { meta, depth -> depth }.collect()
    //PANDEPTH_MERGE
    PANDEPTH_MERGE( ch_all_depth )
    ch_versions = ch_versions.mix(PANDEPTH_MERGE.out.versions)
    // PRODIGAL
    PRODIGAL ( ch_fasta )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    // K-mer
    KMER_COUNT ( ch_fasta, 4  )
    ch_versions = ch_versions.mix(KMER_COUNT.out.versions)

    // Coverage + Kmer + Genes
    SUMMARIZE_FEATURE_MATRIX (
        ch_fasta,
        PANDEPTH_MERGE.out.matrix,
        KMER_COUNT.out.csv,
        PRODIGAL.out.gff
    )

    emit:
    feature_matrix = SUMMARIZE_FEATURE_MATRIX.out.matrix
    versions       = ch_versions
}
