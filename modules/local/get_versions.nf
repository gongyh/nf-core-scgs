process GET_SOFTWARE_VERSIONS {
    output:
    path('software_versions_mqc.yaml'),                           emit: versions

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    trim_galore --version &> v_trim_galore.txt
    bowtie2 --version &> v_bowtie2.txt
    minimap2 -V &> v_minimap2.txt
    samtools --version &> v_samtools.txt
    bedtools --version &> v_bedtools.txt
    preseq &> v_preseq.txt
    qualimap -h &> v_qualimap.txt
    if [ x=`picard MarkDuplicates --version &> v_picard.txt` = 1 ]; then return 0; fi
    #gatk3 -version &> v_gatk.txt
    Rscript -e 'print(packageVersion("AneuFinder"))' &> v_AneuFinder.txt
    spades.py --version &> v_spades.txt
    canu -version &> v_canu.txt
    blastn -version &> v_blast.txt
    quast.py --version &> v_quast.txt
    multiqc --version &> v_multiqc.txt
    diamond version &> v_diamond.txt
    kraken --version | grep Kraken &> v_kraken.txt
    head -n 1 /opt/conda/envs/nf-core-gongyh-scgs/lib/python3.6/site-packages/checkm/VERSION &> v_checkm.txt
    prokka -v &> v_prokka.txt
    emapper.py --version | grep emapper &> v_eggnogmapper.txt
    echo 'v0.0.1' > v_monovar.txt
    blobtools -v &> v_blobtools.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}
