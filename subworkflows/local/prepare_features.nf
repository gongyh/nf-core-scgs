//include { CONTIG_COVERAGE          } from '../../modules/local/pandepth'
include { CONTIG_COVERAGE            } from '../../modules/local/samtools_coverage_combined'
include { PRODIGAL                   } from '../../modules/local/prodigal'
include { KMER_COUNT                 } from '../../modules/local/kmer_count'


workflow PREPARE_FEATURES {
    take:
    ch_fasta
    ch_fai
    ch_bam_for_coverage 

    main:
    ch_versions = Channel.empty()
    /*
    //PANDEPTH_MERGE
    def meta = [id:'merged']
    ch_bam_input = Channel.of( [meta, ch_merged_bam, []] )
    ch_fasta_path = ch_fasta.map { m, file -> file }
    ch_fai_path = ch_fai.map { m, file -> file }
    CONTIG_COVERAGE( ch_bam_input, ch_fasta_path, ch_fai_path )
    ch_depth = CONTIG_COVERAGE.out.depth 
    ch_versions = ch_versions.mix(MERGE_COVERAGE.out.versions)
    ch_coverage_mqc = CONTIG_COVERAGE.out.mqc_tsv
    */
    //samtools
    ch_fasta_path = ch_fasta.map { m, file -> file }
    ch_fai_path = ch_fai.map { m, file -> file }
    CONTIG_COVERAGE( ch_bam_for_coverage, ch_fasta_path, ch_fai_path )

    ch_depth = CONTIG_COVERAGE.out.depth 
    ch_coverage = ch_depth.map { m, depth -> depth }
    ch_coverage_mqc = CONTIG_COVERAGE.out.mqc_tsv
    ch_versions = ch_versions.mix(CONTIG_COVERAGE.out.versions)
    // PRODIGAL
    PRODIGAL ( ch_fasta )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    // K-mer
    KMER_COUNT ( ch_fasta, 4  )
    ch_versions = ch_versions.mix(KMER_COUNT.out.versions)

    emit:
    feature_matrix = Channel.empty()
    coverage_matrix = ch_coverage
    coverage_mqc    = CONTIG_COVERAGE.out.mqc_tsv
    versions        = ch_versions
}
