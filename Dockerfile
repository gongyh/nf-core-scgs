FROM continuumio/miniconda3:4.5.4
LABEL authors="Yanhai Gong" \
      description="Docker image containing all requirements for gongyh/nf-core-scgs pipeline"

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps libpng16-16 cmake gcc g++ hmmer2 && apt-get clean -y 

# Install locale en_US.UTF-8 used by Picard
RUN apt-get install -y locales && sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen && locale-gen

# Install resfinder and pointFinder
RUN cd /opt && git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
RUN cd /opt && git clone https://bitbucket.org/genomicepidemiology/pointfinder.git

# Keep a copy of current pipeline to container
COPY . /opt/nf-core-scgs/

# Install conda environments
COPY environment.yml /
RUN conda install python=3.6 conda=4.6.12 && conda env update -n base -f /environment.yml && conda clean -a

# Install Bioconductor packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('GenomeInfoDbData'); BiocManager::install('AneuFinder')"

# Download silva and busco for Quast 5.x, update taxa db for krona
RUN quast-download-silva && ktUpdateTaxonomy.sh #&& quast-download-busco

# Download GATK and setting
RUN wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-0-ge9d806836" -O GenomeAnalysisTK-3.8.tar.bz2 && \
    tar xjvf GenomeAnalysisTK-3.8.tar.bz2 && gatk3-register GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar && rm -rf ./GenomeAnalysisTK-3.8*

# Install conda environment py27
COPY py27_env.yml /
RUN conda env create -n py27 -f /py27_env.yml && conda clean -a

# Install ACDC
RUN git clone https://github.com/mlux86/acdc.git /tmp/acdc && cd /tmp/acdc && mkdir build && cd build && cmake .. -DBOOST_ROOT=/opt/conda/ && \
    make -j $(nproc) && make install && rm -rf /tmp/acdc && mkdir /acdc

# Install rnammer
RUN mkdir /tmp/rnammer
ADD rnammer-1.2.src.tar.Z /tmp/rnammer
ADD rnammer.patch /tmp/rnammer
RUN cd /tmp/rnammer && tar xf rnammer-1.2.src.tar.Z && rm rnammer-1.2.src.tar.Z && patch < rnammer.patch && mkdir /usr/local/share/rnammer && \
    cp -r * /usr/local/share/rnammer && ln -s /usr/local/share/rnammer/rnammer /usr/local/bin/rnammer && rm -rf /tmp/rnammer

# Build database for blobtools
RUN [ "/bin/bash", "-c", "source activate py27 && blobtools-build_nodesdb && source deactivate" ]

# Install kofamscan
RUN cd /opt && wget ftp://ftp.genome.jp/pub/tools/kofamscan/kofamscan-1.1.0.tar.gz -O kofamscan-1.1.0.tar.gz && \
    tar --no-same-owner -xzvf kofamscan-1.1.0.tar.gz && rm -rf kofamscan-1.1.0.tar.gz

# Clean up
RUN apt-get autoremove --purge && apt-get clean && apt-get autoremove
RUN conda clean -y -a && rm -rf /opt/conda/pkgs/*

# Add default container command
CMD ["/opt/conda/bin/nextflow", "run", "/opt/nf-core-scgs/main.nf", "--help"]

