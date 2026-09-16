nextflow.enable.types = true

//include { CONTIG_COVERAGE          } from '../../modules/local/pandepth'
include { CONTIG_COVERAGE            } from '../../modules/local/samtools_coverage_combined'
include { PRODIGAL                   } from '../../modules/local/prodigal'
include { KMER_COUNT                 } from '../../modules/local/kmer_count'


workflow PREPARE_FEATURES_SINGLE {
    take:
    ch_fasta: Channel<Tuple<Map,Path>>
    ch_fai: Channel<Tuple<Map,Path>>
    ch_bam_for_coverage: Channel<Tuple<Map,Path,List<Path>>>

    main:
    ch_published = channel.empty()
    /*
    //PANDEPTH_MERGE
    def meta = [id:'merged']
    ch_bam_input = channel.of( [meta, ch_merged_bam, []] )
    ch_fasta_path = ch_fasta.map { m, file -> file }
    ch_fai_path = ch_fai.map { m, file -> file }
    CONTIG_COVERAGE( ch_bam_input, ch_fasta_path, ch_fai_path )
    ch_depth = CONTIG_COVERAGE.out.depth
    ch_coverage_mqc = CONTIG_COVERAGE.out.mqc_tsv
    */
    //samtools
    ch_fasta_path = ch_fasta.map { m, file -> file }
    ch_fai_path = ch_fai.map { m, file -> file }
    ch_coverage_input = ch_bam_for_coverage
        .combine(ch_fasta_path)
        .combine(ch_fai_path)
    contig_coverage = CONTIG_COVERAGE(ch_coverage_input)
    ch_coverage = contig_coverage.map { result -> result.depth }
    ch_coverage_mqc = contig_coverage.map { result -> result.mqc_tsv }
    ch_versions = contig_coverage.map { result -> result.versions }
    ch_published = ch_published.mix(contig_coverage.map { result -> [destination: 'coverage_depth', files: result.depth] })
    // PRODIGAL
    prodigal = PRODIGAL(ch_fasta)
    ch_versions = ch_versions.mix(prodigal.map { result -> result.versions })
    ch_published = ch_published.mix(prodigal.map { result -> [destination: 'prodigal', files: result] })

    // K-mer
    kmer_count = KMER_COUNT(ch_fasta, 4)
    ch_versions = ch_versions.mix(kmer_count.map { result -> result.versions })
    ch_published = ch_published.mix(kmer_count.map { result -> [destination: 'kmer', files: result.csv] })

    emit:
    feature_matrix: Channel<Path> = channel.empty()
    coverage_matrix: Channel<Path> = ch_coverage
    coverage_mqc: Channel<Path> = ch_coverage_mqc
    versions: Channel<Path> = ch_versions
    published: Channel<Map> = ch_published
}
