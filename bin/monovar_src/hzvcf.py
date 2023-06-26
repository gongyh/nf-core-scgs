"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Hamim Zafar and Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import time
import os

VCF_META_TEMPLATE = """##fileformat=VCFv4.1
##fileDate={_d.time.tm_year}:{_d.time.tm_mon}:{_d.time.tm_mday}-{_d.time.tm_hour}:{_d.time.tm_min}:{_d.time.tm_sec}
##source=MonoVar_NB
{_d.filter_meta}
{_d.info_meta}
{_d.format_meta}
{_d.ref_meta}
{_d.contig_meta}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{_d.files_meta}
"""

class VCFDocument():

    def __init__(self, out_file, bam_id_list, ref_file):
        self.time = time.localtime()
        self.outf = out_file
        self.outf_rec = open('{}.rec'.format(out_file), 'w')

        self.files_meta = '\t'.join(bam_id_list)
        if ref_file:
            self.ref_meta = '##reference=file:{}'.format(ref_file)
        else:
            self.ref_meta = '##reference=?'
        
        info_fields = [
            ('AC', 'A', 'Integer',
                'Allele count in genotypes, for each ALT allele, in the same order as listed'),
            ('AF', 'A', 'Float',
                'Allele Frequency, for each ALT allele, in the same order as listed'),
            ('AN', '1', 'Integer',
                'Total number of alleles in called genotypes'),
            ('BaseQRankSum', '1', 'Float',
                'Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities'),
            ('DP', '1', 'Integer',
                'Approximate read depth; some reads may have been filtered'),
            ('QD', '1', 'Float',
                'Variant Confidence/Quality by Depth'),
            ('SOR', '1', 'Float',
                'Symmetric Odds Ratio of 2x2 contingency table to detect strand bias'),
            ('MPR', '1', 'Float',
                'Log Odds Ratio of maximum value of probability of observing non-ref allele to the probability of observing zero non-ref allele'),
            ('PSARR', '1', 'Float',
                'Ratio of per-sample Alt allele supporting reads to Ref allele supporting reads')
        ]
        self.info_meta = '\n'.join([
            '##INFO=<ID={},Number={},Type={},Description="{}">'.format(*i) \
                for i in info_fields])

        filter_fields = [
            ('LowQual', 'Low quality'),
            ('NoConsensus', 'Only 1 sample contains called SNV')
        ]
        self.filter_meta = '\n'.join([
            '##FILTER=<ID={},Description="{}">'.format(*i) \
                for i in filter_fields])

        format_fields = [
            ('AD', '.', 'Integer', 
                'Allelic depths for the ref and alt alleles in the order listed'),
            ('DP', '1', 'Integer',
                'Approximate read depth (reads with MQ=255 or with bad mates are filtered)'),
            ('GQ', '1', 'Integer', 'Genotype Quality'),
            ('GT', '1', 'String', 'Genotype'),
            ('PL', 'G', 'Integer',
                'Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification')
        ]
        self.format_meta = '\n'.join([
            '##FORMAT=<ID={},Number={},Type={},Description="{}">'.format(*i) \
                for i in format_fields])

        self.contig_meta = '##contig=<ID=-1,eta=-1>'


    def append_record(self, rec_data):
        rec_str = '\t'.join(rec_data) + '\n'
        self.outf_rec.write(rec_str)


    def add_contigs(self, contigs):
        self.contig_meta = '\n'.join(['##contig=<ID={},eta=-1>'.format(i) \
            for i in contigs])


    def add_header(self):
        with open(self.outf, 'w') as f_out:
            f_out.write(VCF_META_TEMPLATE.format(_d=self))
            with open(self.outf_rec.name, 'r') as f_rec:
                f_out.write(f_rec.read())

        os.remove(self.outf_rec.name)


    def close_records(self):
        self.outf_rec.close()
        

if __name__ == '__main__':
    print('Here be dragons...')
