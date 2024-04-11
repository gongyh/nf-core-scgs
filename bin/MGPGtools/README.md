# MGPGtools
Tools for microbial graph pangenomes

## Contents

## Installation
### Pre-installation
#### Required Python packages
+ python (>=3.4)
+ prettytable
+ gfapy
+ Biopython
+ toytree
+ pandas
+ numpy

#### External requirements
+ odgi (0.8.3)
+ seqkit (v2.5.1)
+ clustalW (v2.1)
+ aster (v1.16)
+ gffread (v0.12.7)

### Quick environment setup
```
conda/mamba/micromamba env create -f environment.yml
conda/mamba/micromamba activate mgpg-tools
python main.py -h
```

## Arguments

## Examples
    usage: python main.py stat -db DATABASE -rank [domain|phylum|class|order|family|genus] -name NAME -outdir OUTDIR [-h Help]
           eg. python main.py stat -db path_to_database -rank order -name Lactobacillales -outdir path_to_dir
           output: Return a TSV file containing all metagenomes from the pan-genome database with a specific rank and specific genus.
          |  Domain  |   Phylum  |  Class  |      Order      |       Family      |     Genus    | Nodes | Links | genomes |       Ref       |
          | ---------| ----------|---------|-----------------|-------------------|--------------|-------|-------|---------|-----------------|
          | Bacteria | Bacillota | Bacilli | Lactobacillales | Carnobacteriaceae | Alkalibacter |  5915 |  7643 |    16   | GCA_003446375.1 |

           python main.py search -db DATABASE -name NAME [-gene|-cds] [-h]
           eg. python main.py search -db path_to_database -name Listeria [-gene LSE_RS00145] [-cds WP_003721652.1]
           output: Return a table containing basic information of metagenomes from a specific genus in the pan-genome database, along with information on specific genes/CDS.
          |   name   |       ref       | chromosomes |   size  |            path             | start |  end  |
          |----------|-----------------|-------------|---------|-----------------------------|-------|-------|
          | Listeria | GCF_000027145.1 | NC_013891.1 | 2797636 | GCF_000027145#1#NC_013891.1 | 31414 | 31752 |

           python main.py viz -db DATABASE -name NAME -outdir OUTDIR [-outName OUTNAME] [-t Thread] [-x WIDTH] [-y HEIGTH] [-r PATH_RANGE] [-h Help]
           eg. python main.py viz -db path_to_database -name Listeria -outdir path_to_dir -outName Listeria -t 16 -x 300 -y 200 [-r GCF_000027145#1#NC_013891.1:31414-31752]
           output: Return a PNG image containing the graphical representation of the pan-genome at a specific location.

           python main.py describe -db DATABASE -name NAME -outdir OUTDIR [-t Thread] [-h Help]
           eg. python main.py describe -db path_to_database -name Listeria -outdir path_to_dir -t 16
           output: Return a TSV file containing information on the pan-genome.

           python main.py tree [-gfa GFA] [-sampleTxt sampleTxt] [-label LABEL] -db DATABASE -name NAME -gene GENE -outdir OUTDIR [-t Thread] [-h Help]
           eg. python main.py tree -db path_to_database -name Listeria -gene LSE_RS00145 -outdir path_to_dir -t 16
           eg. python main.py tree -gfa test/sample.gfa -sampleTxt test/sample.txt -label sample -db path_to_database -name Listeria -gene LSE_RS00145 -outdir path_to_dir -t 16
           output: Return an SVG format vector graphic that includes the phylogenetic tree depicting the evolution of the gene in the graph pangenome.

           python main.py core [-gfa GFA] -db DATABASE -name NAME -outdir OUTDIR [-t Thread] [-h Help]
           eg. python main.py core -db path_to_database -name Moorella -outdir path_to_dir -t 16
           output: Return four TXT file(0-15%.txt, 15-95%.txt, 95-99%.txt, 99-100%.txt, 100%.txt), Gene names with various threshold values included.
