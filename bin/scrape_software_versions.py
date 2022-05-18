#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
from collections import OrderedDict
import re

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    'gongyh/nf-core-scgs': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'Trim Galore!': ['v_trim_galore.txt', r"version (\S+)"],
    'Bowtie2': ['v_bowtie2.txt', r"bowtie2-align-s version (\S+)"],
    'Minimap2': ['v_minimap2.txt', r"(\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'BEDTools': ['v_bedtools.txt', r"bedtools v(\S+)"],
    'Preseq': ['v_preseq.txt', r"Version: (\S+)"],
    'QualiMap': ['v_qualimap.txt', r"QualiMap v.(\S+)"],
    'Picard': ['v_picard.txt', r"(\S+)"],
    'Monovar': ['v_monovar.txt', r"v(\S+)"],
    'AneuFinder': ['v_AneuFinder.txt', r'\[1\] ‘(\S+)’'],
    'Spades': ['v_spades.txt', r"SPAdes v(\S+)"],
    'Canu': ['v_canu.txt', r"canu (\S+)"],
    'BLAST': ['v_blast.txt', r"blastn: (\S+)"],
    'Diamond': ['v_diamond.txt', r"diamond version (\S+)"],
    'Kraken': ['v_kraken.txt', r"Kraken version (\S+)"],
    'CheckM': ['v_checkm.txt', r"(\S+)"],
    'Prokka': ['v_prokka.txt', r"prokka (\S+)"],
    'eggNOG-mapper': ['v_eggnogmapper.txt', r"emapper-(\S+)"],
    'BlobTools': ['v_blobtools.txt', r"blobtools v(\S+)"],
    'Quast': ['v_quast.txt', r"QUAST v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
}
results = OrderedDict()
results['gongyh/nf-core-scgs'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['Trim Galore!'] = '<span style="color:#999999;\">N/A</span>'
results['Bowtie2'] = '<span style="color:#999999;\">N/A</span>'
results['Minimap2'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['BEDTools'] = '<span style="color:#999999;\">N/A</span>'
results['Preseq'] = '<span style="color:#999999;\">N/A</span>'
results['QualiMap'] = '<span style="color:#999999;\">N/A</span>'
results['Picard'] = '<span style="color:#999999;\">N/A</span>'
results['Monovar'] = '<span style="color:#999999;\">N/A</span>'
results['AneuFinder'] = '<span style="color:#999999;\">N/A</span>'
results['Spades'] = '<span style="color:#999999;\">N/A</span>'
results['Canu'] = '<span style="color:#999999;\">N/A</span>'
results['BLAST'] = '<span style="color:#999999;\">N/A</span>'
results['Diamond'] = '<span style="color:#999999;\">N/A</span>'
results['Kraken'] = '<span style="color:#999999;\">N/A</span>'
results['CheckM'] = '<span style="color:#999999;\">N/A</span>'
results['Prokka'] = '<span style="color:#999999;\">N/A</span>'
results['eggNOG-mapper'] = '<span style="color:#999999;\">N/A</span>'
results['BlobTools'] = '<span style="color:#999999;\">N/A</span>'
results['Quast'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'gongyh/nf-core-scgs-software-versions'
section_name: 'gongyh/nf-core-scgs Software Versions'
section_href: 'https://github.com/gongyh/nf-core-scgs'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
