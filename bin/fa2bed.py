#!/usr/bin/env python2
from __future__ import print_function

from Bio.Seq import Seq
import sys
from Bio.SeqUtils import GC, GC_skew
from Bio import SeqIO
import click

@click.command()
@click.option('--window', default=100, help='length of sliding window (default: 100bp)')
@click.argument('fa', type=click.File('rb'))
def fa2bed(fa, window):

    # genome bed file
    genomeBed = open("genome.bed","w")
    #genomeBed.write("chr\tstart\tend\n")

    record_dict = SeqIO.to_dict(SeqIO.parse(fa, "fasta"))

    # GC% distribution file
    gcBed = open("gc.bed","w")
    #gcBed.write("chr\tstart\tend\tgc\n")

    #GC-skew bed file
    gcSkewBed = open("gcSkew.bed","w")
    #gcSkewBed.write("chr\tstart\tend\tgcSkew\n")

    window_len = int(window)

    for k,v in record_dict.items():
        chrom = v.id
        chrom_seq = v.seq
        length = len(chrom_seq)
        genomeBed.write(chrom+"\t0\t%d\n"%length)

        for i in range(0, length, window_len):
            start = i
            end = i + window_len
            s = chrom_seq[start : end]
            gc = GC(s)
            skew = GC_skew(s, window_len)[0]
            gcBed.write(chrom+"\t%d\t%d\t%.3f\n"%(start,end,gc))
            gcSkewBed.write(chrom+"\t%d\t%d\t%.3f\n"%(start,end,skew))

    genomeBed.close()
    gcBed.close()
    gcSkewBed.close()

if __name__ == "__main__":
    fa2bed()

