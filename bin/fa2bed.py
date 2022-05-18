#!/usr/bin/env python
from __future__ import print_function

from Bio.Seq import Seq
import sys
from Bio.SeqUtils import GC, GC_skew
from Bio import SeqIO
import click

@click.command()
@click.option('--window', default=10000, help='length of sliding window (default: 10Kbp)')
@click.option('--step', default=200, help='length of step window (default: 200bp)')
@click.argument('fa', type=click.Path(exists=True))
def fa2bed(fa, window, step):

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
    step_len = int(step)

    for k,v in record_dict.items():
        chrom = v.id
        chrom_seq = v.seq
        length = len(chrom_seq)
        genomeBed.write(chrom+"\t0\t%d\n"%length)

        for i in range(0, length, step_len):
            start0 = i
            end0 = i+step_len if i+step_len<=length else length
            start = int(i-window_len/2) if i-window_len/2>=0 else 0
            end = int(i+window_len/2) if i+window_len/2<=length else length
            s = chrom_seq[start : end]
            gc = GC(s)
            skew = GC_skew(s, window_len)[0]
            gcBed.write(chrom+"\t%d\t%d\t%.3f\n"%(start0,end0,gc))
            gcSkewBed.write(chrom+"\t%d\t%d\t%.3f\n"%(start0,end0,skew))

    genomeBed.close()
    gcBed.close()
    gcSkewBed.close()

if __name__ == "__main__":
    fa2bed()

