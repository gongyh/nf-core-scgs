#!/bin/bash

## prepare databases for the pipeline

echo """

The following databases are supported by gongyh/scgs pipeline:
  1) NCBI nt database
  2) Uniprot proteomes database
  3) Kraken database
  4) CheckM database
  5) EggNOG v5 database
  6) KOfam HMM database
  7) ResFinder database
  8) PointFinder database
  9) Funannotate database

Current working directory is ${PWD}.
All database files will be downloaded into ${PWD}.
Please make sure you have ownership for ${PWD}.

"""

read -p "Please choose the database you need: " choice

if [ $choice -eq 1 ]; then
  mkdir -p $PWD/nt
  cd $PWD/nt
  update_blastdb.pl nt
fi

if [ $choice -eq 2 ]; then
  mkdir -p $PWD/uniprot
  cd uniprot
  # Download database
  wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Reference_Proteomes_2017_07.tar.gz
  # Unpack protein FASTAs for each kingdom
  parallel -j4 'gunzip {}' ::: `ls | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional'`
  # Concatenate all protein sequences into uniprot_ref_proteomes.fasta
  cat */*.fasta > uniprot_ref_proteomes.fasta
  # Simplify sequence IDs
  cat uniprot_ref_proteomes.fasta | sed -r 's/(^>sp\|)|(^>tr\|)/>/g' | cut -f1 -d"|" > temp; mv temp uniprot_ref_proteomes.fasta
  # Make Diamond database
  diamond makedb --in uniprot_ref_proteomes.fasta -d uniprot_ref_proteomes.diamond
  # Subset mapping file to only contain NCBI TaxID entries
  cat */*.idmapping | grep "NCBI_TaxID" > uniprot_ref_proteomes.taxids
fi

if [ $choice -eq 3 ]; then
  mkdir -p $PWD/kraken
  cd kraken
  wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171101_8GB_dustmasked.tgz
  tar xzvf minikraken_20171101_8GB_dustmasked.tgz
fi

if [ $choice -eq 4 ]; then
  mkdir -p $PWD/checkm
  cd $PWD/checkm
  wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
  tar -xzvf checkm_data_2015_01_16.tar.gz
fi

if [ $choice -eq 5 ]; then
  mkdir -p $PWD/eggnog
  cd $PWD/eggnog
  conda activate py27
  download_eggnog_data.py -y
  conda deactivate
fi

if [ $choice -eq 6 ]; then
  mkdir -p $PWD/kofam
  cd $PWD/kofam
  wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
  tar xzvf profiles.tar.gz
  wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
  gzip -d ko_list.gz
fi

if [ $choice -eq 7 ]; then
  mkdir -p $PWD/resfinder
  cd $PWD/resfinder
  git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
fi

if [ $choice -eq 8 ]; then
  mkdir -p $PWD/pointFinder
  cd $PWD/pointFinder
  git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git
fi

if [ $choice -eq 9 ]; then
  mkdir -p $PWD/funannotate
  funannotate setup -d $PWD/funannotate
fi


echo "Done!"

