nextflow.enable.types = true

//include { CONTIG_COVERAGE ;MERGE_COVERAGE } from '../../modules/local/pandepth'
include { CONTIG_COVERAGE ;MERGE_COVERAGE } from '../../modules/local/samtools_coverage_combined'
include { PRODIGAL                   } from '../../modules/local/prodigal'
include { KMER_COUNT                 } from '../../modules/local/kmer_count'
include { SUMMARIZE_FEATURE_MATRIX   } from '../../modules/local/summarize_feature_matrix'

workflow PREPARE_FEATURES_MULTI {
    take:
    ch_fasta: Channel<Tuple<Map,Path>>
    ch_fai: Channel<Tuple<Map,Path>>
    ch_bams: Channel<Tuple<Map,Path>>

    main:
    ch_versions = channel.empty()
    /*
    //PANDEPTH_MERGE
    ch_fasta_file = ch_fasta
        .map { item -> item instanceof List ? item : [item] }
        .flatten()
        .filter { path -> path.toString().endsWith('.fasta') || path.toString().endsWith('.fa') }
        .first()
    ch_fai_file = ch_fai
        .map { item -> item instanceof List ? item : [item] }
        .flatten()
        .filter { path -> path.toString().endsWith('.fai') }
        .first()
    ch_bams_with_bai = ch_bams.map { meta, bam -> [meta, bam, []] }
    ch_depth = CONTIG_COVERAGE( ch_bams_with_bai, ch_fasta_file, ch_fai_file ).depth
    ch_all_depth = ch_depth.map { meta, depth -> depth }.collect()

    MERGE_COVERAGE( ch_all_depth )
    ch_versions = ch_versions.mix(MERGE_COVERAGE.out.versions)
    */
    //samtools
    ch_fasta_file = ch_fasta.map { _meta, fasta -> fasta }
    ch_fai_file = ch_fai.map { _meta, fai -> fai }

    ch_bams_with_bai = ch_bams.map { meta, bam -> tuple(meta, bam, [] as List<Path>) }

    ch_coverage_input = ch_bams_with_bai
        .combine(ch_fasta_file)
        .combine(ch_fai_file)
    contig_coverage = CONTIG_COVERAGE(ch_coverage_input)
    ch_all_depth = contig_coverage.map { result -> result.depth }.collect()

    merge_coverage = MERGE_COVERAGE(ch_all_depth)
    ch_coverage_mqc = contig_coverage.map { result -> result.mqc_tsv }
    ch_versions = ch_versions.mix(merge_coverage.map { result -> result.versions })
    // PRODIGAL
    prodigal = PRODIGAL(ch_fasta)
    ch_versions = ch_versions.mix(prodigal.map { result -> result.versions })

    // K-mer
    kmer_count = KMER_COUNT(ch_fasta, 4)
    ch_versions = ch_versions.mix(kmer_count.map { result -> result.versions })

    // Coverage + Kmer + Genes
    ch_feature_input = ch_fasta
        .combine(merge_coverage.map { result -> result.matrix })
        .combine(kmer_count.map { result -> tuple(result.kmer, result.csv) })
        .combine(prodigal.map { result -> result.gff })
    summarize_feature_matrix = SUMMARIZE_FEATURE_MATRIX(ch_feature_input)

    emit:
    feature_matrix: Channel<Path> = summarize_feature_matrix.map { result -> result.matrix }
    coverage_matrix: Value<Path> = merge_coverage.map { result -> result.matrix }
    coverage_mqc: Channel<Path> = ch_coverage_mqc
    versions: Channel<Path> = ch_versions
}
