#!/bin/bash

## prepare databases for the pipeline

echo """

The following databases are supported by gongyh/scgs pipeline:
  1) NCBI nt database (latest version)
  2) Uniprot proteomes database (Reference_Proteomes_2017_07)
  3) Kraken database (minikraken_20171101_8GB_dustmasked)
  4) CheckM database (checkm_data_2015_01_16, no longer needed after v1.0-rc)
  5) EggNOG v5 database (v5)
  6) KOfam HMM database (latest)
  7) ResFinder database (latest)
  8) PointFinder database (latest)
  9) Funannotate database (latest, only needed for eukaryotic microbes, e.g. fungi)
  10) PlasmidFinder database (latest)
  11) VirulenceFinder database (latest)
  12) Pfam-A database (Pfam31.0)
  13) EukCC database (V1.1, 20191023_1)

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
  . /opt/conda/etc/profile.d/conda.sh
  conda activate scgs_py27
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
  . /opt/conda/etc/profile.d/conda.sh
  conda activate scgs_py27
  mkdir -p $PWD/funannotate
  funannotate setup -d $PWD/funannotate
  conda deactivate
fi

if [ $choice -eq 10 ]; then
  mkdir -p $PWD/plasmidfinder
  cd $PWD/plasmidfinder
  git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git
fi

if [ $choice -eq 11 ]; then
  mkdir -p $PWD/virulencefinder
  cd $PWD/virulencefinder
  git clone https://bitbucket.org/genomicepidemiology/virulencefinder_db.git
fi

if [ $choice -eq 12 ]; then
  mkdir -p $PWD/PFAMDB
  cd $PWD/PFAMDB
  curl ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz > Pfam-A.hmm.gz
  gzip -d Pfam-A.hmm.gz
  hmmpress Pfam-A.hmm 1>/dev/null 2>/dev/null
  echo "Please set system environment variable PFAMDB_PATH to $PWD/PFAMDB before using ezTree."
fi

if [ $choice -eq 13 ]; then
  wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc_db_v1.1.tar.gz
  tar -xzvf eukcc_db_v1.1.tar.gz
  mv eukcc_db_20191023_1 eukcc_db
fi

echo "Done!"

