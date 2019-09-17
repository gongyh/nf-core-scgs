Bootstrap: docker
FROM: continuumio/miniconda3:4.5.4
IncludeCmd: yes

%help
Singularity container for gongyh/scgs pipeline.

%files
environment.yml /
py27_env.yml /

%labels
Author gongyh@qibebt.ac.cn

%post
export PATH="/opt/conda/bin:$PATH"
apt-get update && apt-get install -y procps libpng16-16 && apt-get clean -y
apt-get install -y locales && sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen && locale-gen
conda install python=3.6 conda=4.6.12 && conda env update -n base -f /environment.yml && conda clean -a
R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('GenomeInfoDbData'); BiocManager::install('AneuFinder')"
quast-download-silva && ktUpdateTaxonomy.sh
wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-0-ge9d806836" -O GenomeAnalysisTK-3.8.tar.bz2 && \
    tar xjvf GenomeAnalysisTK-3.8.tar.bz2 && gatk3-register GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar && rm -rf ./GenomeAnalysisTK-3.8*
conda env create -n py27 -f /py27_env.yml && conda clean -a
/bin/bash -c "source activate py27 && blobtools-build_nodesdb && source deactivate"
cd /opt && git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
cd /opt && git clone https://bitbucket.org/genomicepidemiology/pointfinder.git
apt-get autoremove --purge && apt-get clean && apt-get autoremove
conda clean -y -a && rm -rf /opt/conda/pkgs/*
