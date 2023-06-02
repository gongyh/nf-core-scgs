process POINTFINDER {
    tag "$prefix"
    publishDir "${params.outdir}/ARG", mode: 'copy'

    input:
    path contigs
    path db

    output:
    path("${prefix}/*")

    when:
    !euk && params.point && pointfinder_db

    script:
    prefix = contigs.toString() - ~/(\.ctgs\.fasta)?(\.ctgs)?(\.fasta)?(\.fa)?$/
    def species = params.pointfinder_species
    def known_snp = params.only_known ? "" : "-l 0.4 -r all -u"
    """
    mkdir -p $prefix
    python /opt/pointfinder/PointFinder.py -p $db \
    -m blastn -m_p /opt/conda/bin/blastn $known_snp \
    -i $contigs -o $prefix -s $species
    rm -rf $prefix/tmp
    """
}
