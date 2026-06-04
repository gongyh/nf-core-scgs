include { TNF_RPKM } from '../../modules/local/dcvbin/extract_tnf_rpkm_feature/main'
include { CONTIG_EMBEDDING } from '../../modules/local/dcvbin/extract_fpf_feature/main'
include { FEATURE_FUSION } from '../../modules/local/dcvbin/feature_fusion_vae/main'
include { CONTIG_KMER } from '../../modules/local/dcvbin/extract_kmer_feature/main'
include { MARKER_NCLUSTERS } from '../../modules/local/dcvbin/nclusters_marker/main'
include { DCVBIN_BIN } from '../../modules/local/dcvbin/clustering_bins/main'

workflow DCVBIN {
    take:
    contigs_fasta // contigs >= 2k
    sorted_bam
    
    main:
    
    // contig features from DNABERT-S
    CONTIG_EMBEDDING(contigs_fasta, params.DNABERTS_dir)
    
    // TNF & RPKM
    TNF_RPKM(contigs_fasta, sorted_bam)
    
    // VAE feature fusion
    FEATURE_FUSION(
        CONTIG_EMBEDDING.out.fpf,
        TNF_RPKM.out.tnf,
        TNF_RPKM.out.rpkm
    )
    
    // k-mer feature
    CONTIG_KMER(contigs_fasta)
    
    // Initial number of clusters by marker genes
    MARKER_NCLUSTERS(CONTIG_KMER.out.kmer, contigs_fasta)
    
    // binning
    DCVBIN_BIN(
        FEATURE_FUSION.out.features,
        MARKER_NCLUSTERS.out.marker_cv,
        contigs_fasta
    )
    
    emit:
    bins_dir = DCVBIN_BIN.out.bins_dir
    label_file = DCVBIN_BIN.out.label_file
}
