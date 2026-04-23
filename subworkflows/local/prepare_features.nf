include { SAMTOOLS_COVERAGE_COMBINED } from '../../modules/local/samtools_coverage_combined'
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

    all_bams = ch_bams.map { it[1] }.collect()
    fasta = ch_fasta.map { meta, fasta_file -> fasta_file }.first()
    fai = ch_fai
        .map { it instanceof List ? it.flatten() : [it] }
        .flatten()
        .filter { it.toString().endsWith('.fai') }
        .first()
    // COVERAGE_COMBINED
    SAMTOOLS_COVERAGE_COMBINED( all_bams, fasta, fai )
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE_COMBINED.out.versions)
    // PRODIGAL
    PRODIGAL ( ch_fasta )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    // K-mer
    KMER_COUNT ( ch_fasta, 4  )
    ch_versions = ch_versions.mix(KMER_COUNT.out.versions)

    // Coverage + Kmer + Genes
    SUMMARIZE_FEATURE_MATRIX (
        ch_fasta,
        SAMTOOLS_COVERAGE_COMBINED.out.matrix,
        KMER_COUNT.out.csv,
        PRODIGAL.out.gff
    )

    emit:
    feature_matrix = SUMMARIZE_FEATURE_MATRIX.out.matrix
    versions       = ch_versions
}
