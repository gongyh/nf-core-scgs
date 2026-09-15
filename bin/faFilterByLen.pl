#!/usr/bin/env perl
use warnings;
use strict;

use Bio::SeqIO;

my $contig = $ARGV[0];
my $minLen = $ARGV[1];

my $ctgs  = Bio::SeqIO->new( -format => 'fasta', -file => $contig);

while( my $ctg = $ctgs->next_seq() ) {
    my $id  = $ctg->primary_id; chomp $id;
    my $seq = $ctg->seq; chomp $seq;
    my $len = length($seq);
    if( $len >= $minLen ) {
        print ">$id","\n",$seq,"\n";
    }
}
