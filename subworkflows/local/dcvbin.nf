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
    ch_published = channel.empty()
    // contig features from DNABERT-S
    ch_contigs_fasta = contigs_fasta.map { meta, fasta ->
        tuple(meta, fasta)
    }
    ch_embedding_input = ch_contigs_fasta.map { meta, fasta ->
        tuple(meta, fasta, file(params.DNABERTS_dir, type: 'dir'))
    }
    contig_embedding = CONTIG_EMBEDDING(ch_embedding_input)
    ch_published = ch_published.mix(contig_embedding.map { result -> [destination: 'dcvbin_embeddings', files: result] })

    // TNF & RPKM
    ch_tnf_rpkm_input = ch_contigs_fasta.combine(sorted_bam)
    tnf_rpkm = TNF_RPKM(ch_tnf_rpkm_input)
    ch_published = ch_published.mix(tnf_rpkm.map { result -> [destination: "dcvbin_tnf_rpkm/${result.meta.id}", files: result] })

    // VAE feature fusion
    ch_feature_fusion_input = contig_embedding
        .combine(tnf_rpkm)
        .map { embedding, tnf_rpkm_result ->
            tuple(embedding.meta, embedding.fpf, tnf_rpkm_result.tnf, tnf_rpkm_result.rpkm)
        }
    feature_fusion = FEATURE_FUSION(ch_feature_fusion_input)
    ch_published = ch_published.mix(feature_fusion.map { result -> [destination: "dcvbin_vae/${result.meta.id}", files: result] })

    // k-mer feature
    contig_kmer = CONTIG_KMER(ch_contigs_fasta)
    ch_published = ch_published.mix(contig_kmer.map { result -> [destination: 'dcvbin_kmer', files: result] })

    // Initial number of clusters by marker genes
    ch_marker_input = contig_kmer
        .map { result -> tuple(result.meta, result.kmer) }
        .combine(ch_contigs_fasta.map { _meta, fasta -> fasta })
    marker_nclusters = MARKER_NCLUSTERS(ch_marker_input)
    ch_published = ch_published.mix(marker_nclusters.map { result -> [destination: 'dcvbin_marker', files: result] })

    // binning
    ch_binning_input = feature_fusion
        .map { result -> tuple(result.meta, result.features) }
        .combine(marker_nclusters.map { result -> result.marker_cv })
        .combine(ch_contigs_fasta.map { _meta, fasta -> fasta })
    dcvbin_bin = DCVBIN_BIN(ch_binning_input)
    ch_published = ch_published.mix(dcvbin_bin.map { result -> [destination: 'dcvbin_bins', files: result] })

    emit:
    bins_dir: Channel<Tuple<Map,Path>> = dcvbin_bin.map { result -> tuple(result.meta, result.bins_dir) }
    label_file: Channel<Tuple<Map,Path>> = dcvbin_bin.map { result -> tuple(result.meta, result.label_file) }
    scaffolds2bin: Channel<Tuple<Map,Path>> = dcvbin_bin.map { result -> tuple(result.meta, result.scaffolds2bin) }
    mqc_tsv: Channel<Tuple<Map,Path>> = dcvbin_bin.map { result -> tuple(result.meta, result.mqc_tsv) }
    published: Channel<Map> = ch_published
}
