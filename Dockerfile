FROM continuumio/miniconda3:4.5.4
LABEL authors="Yanhai Gong" \
      description="Docker image containing all requirements for gongyh/nf-core-scgs pipeline"

COPY environment.yml /
RUN conda env update -n base -f /environment.yml && conda clean -a

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y 

# Install Bioconductor packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('GenomeInfoDbData'); BiocManager::install('AneuFinder')"

# py27
RUN conda create -n py27 -c bioconda -y python=2.7 checkm-genome biopython click monovar && conda clean -a
RUN mkdir -p /opt/checkm-data && cd /opt/checkm-data && \
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar xzvf checkm_data_2015_01_16.tar.gz && rm -rf checkm_data_2015_01_16.tar.gz
RUN [ "/bin/bash", "-c", "source activate py27 && ((echo /opt/checkm-data; sleep 2; echo /opt/checkm-data) | checkm data setRoot) && source deactivate" ]

# clean up
RUN apt-get autoremove --purge && apt-get clean && apt-get autoremove
RUN conda clean -y -a && rm -rf /opt/conda/pkgs/*
