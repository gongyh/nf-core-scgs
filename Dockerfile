FROM continuumio/miniconda3:4.5.4
LABEL authors="Yanhai Gong" \
      description="Docker image containing all requirements for gongyh/nf-core-scgs pipeline"

COPY environment.yml /
RUN conda env update -n base -f /environment.yml && conda clean -a

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y 

# Install locale en_US.UTF-8 used by Picard
RUN apt-get install -y locales && sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen && locale-gen

# Install Bioconductor packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('GenomeInfoDbData'); BiocManager::install('AneuFinder')"

# Fix circos gd problem
#RUN cd /opt/conda/lib && ln -s libwebp.so.6 libwebp.so.7

# Download silva and busco for Quast 5.x
RUN quast-download-silva && quast-download-busco

# Download GATK and setting
RUN wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-0-ge9d806836" -O GenomeAnalysisTK-3.8.tar.bz2 && \
    tar xjvf GenomeAnalysisTK-3.8.tar.bz2 && gatk3-register GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar && rm -rf ./GenomeAnalysisTK-3.8*

# py27
RUN conda create -n py27 -c bioconda -y python=2.7 checkm-genome biopython click monovar blobtools=1.0.1 && conda clean -a
RUN mkdir -p /opt/checkm-data && cd /opt/checkm-data && \
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar xzvf checkm_data_2015_01_16.tar.gz && rm -rf checkm_data_2015_01_16.tar.gz
RUN [ "/bin/bash", "-c", "source activate py27 && ((echo /opt/checkm-data; sleep 2; echo /opt/checkm-data) | checkm data setRoot) && source deactivate" ]
RUN [ "/bin/bash", "-c", "source activate py27 && blobtools-build_nodesdb && source deactivate" ]

# clean up
RUN apt-get autoremove --purge && apt-get clean && apt-get autoremove
RUN conda clean -y -a && rm -rf /opt/conda/pkgs/*
