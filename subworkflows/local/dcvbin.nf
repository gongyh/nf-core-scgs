nextflow.enable.types = true

include { TNF_RPKM } from '../../modules/local/dcvbin/extract_tnf_rpkm_feature/main'
include { CONTIG_EMBEDDING } from '../../modules/local/dcvbin/extract_fpf_feature/main'
include { FEATURE_FUSION } from '../../modules/local/dcvbin/feature_fusion_vae/main'
include { CONTIG_KMER } from '../../modules/local/dcvbin/extract_kmer_feature/main'
include { MARKER_NCLUSTERS } from '../../modules/local/dcvbin/nclusters_marker/main'
include { DCVBIN_BIN } from '../../modules/local/dcvbin/clustering_bins/main'

workflow DCVBIN {
    take:
    contigs_fasta: Channel<Tuple<Map,Path>> // contigs >= 2k
    sorted_bam: Channel<Path>

    main:
    // contig features from DNABERT-S
    ch_embedding_input = contigs_fasta.map { meta, fasta ->
        tuple(meta, fasta, file(params.DNABERTS_dir, type: 'dir'))
    }
    contig_embedding = CONTIG_EMBEDDING(ch_embedding_input)

    // TNF & RPKM
    ch_tnf_rpkm_input = contigs_fasta.combine(sorted_bam)
    tnf_rpkm = TNF_RPKM(ch_tnf_rpkm_input)

    // VAE feature fusion
    ch_feature_fusion_input = contig_embedding
        .combine(tnf_rpkm)
        .map { embedding, tnf_rpkm_result ->
            tuple(embedding.meta, embedding.fpf, tnf_rpkm_result.tnf, tnf_rpkm_result.rpkm)
        }
    feature_fusion = FEATURE_FUSION(ch_feature_fusion_input)

    // k-mer feature
    contig_kmer = CONTIG_KMER(contigs_fasta)

    // Initial number of clusters by marker genes
    ch_marker_input = contig_kmer
        .map { result -> tuple(result.meta, result.kmer) }
        .combine(contigs_fasta.map { _meta, fasta -> fasta })
    marker_nclusters = MARKER_NCLUSTERS(ch_marker_input)

    // binning
    ch_binning_input = feature_fusion
        .map { result -> tuple(result.meta, result.features) }
        .combine(marker_nclusters.map { result -> result.marker_cv })
        .combine(contigs_fasta.map { _meta, fasta -> fasta })
    dcvbin_bin = DCVBIN_BIN(ch_binning_input)

    emit:
    bins_dir: Channel<Tuple<Map,Path>> = dcvbin_bin.map { result -> tuple(result.meta, result.bins_dir) }
    label_file: Channel<Tuple<Map,Path>> = dcvbin_bin.map { result -> tuple(result.meta, result.label_file) }
    scaffolds2bin: Channel<Tuple<Map,Path>> = dcvbin_bin.map { result -> tuple(result.meta, result.scaffolds2bin) }
    mqc_tsv: Channel<Tuple<Map,Path>> = dcvbin_bin.map { result -> tuple(result.meta, result.mqc_tsv) }
    versions: Channel<Path> = dcvbin_bin.map { result -> result.versions }
}
