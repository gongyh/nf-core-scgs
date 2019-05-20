Bootstrap: docker
FROM: continuumio/miniconda3:4.5.4
IncludeCmd: yes

%help
Singularity container for gongyh/scgs pipeline.

%files
    environment.yml /

%labels
    Author gongyh@qibebt.ac.cn

%post
    apt-get update && apt-get install -y procps && apt-get clean -y
    apt-get install -y locales && sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen && locale-gen
    conda install conda=4.6.12 && conda env update -n base -f /environment.yml && conda clean -a
    R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('GenomeInfoDbData'); BiocManager::install('AneuFinder')"
    quast-download-silva && quast-download-busco
    wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-0-ge9d806836" -O GenomeAnalysisTK-3.8.tar.bz2 && \
      tar xjvf GenomeAnalysisTK-3.8.tar.bz2 && gatk3-register GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar && rm -rf ./GenomeAnalysisTK-3.8*
    conda create -n py27 -c bioconda -y python=2.7 checkm-genome biopython click monovar blobtools=1.0.1 && conda clean -a
    mkdir -p /opt/checkm-data && cd /opt/checkm-data && \
      wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
      tar xzvf checkm_data_2015_01_16.tar.gz && rm -rf checkm_data_2015_01_16.tar.gz
    /bin/bash -c "source activate py27 && ((echo /opt/checkm-data; sleep 2; echo /opt/checkm-data) | checkm data setRoot) && source deactivate"
    /bin/bash -c "source activate py27 && blobtools-build_nodesdb && source deactivate"
    apt-get autoremove --purge && apt-get clean && apt-get autoremove
    conda clean -y -a && rm -rf /opt/conda/pkgs/*
