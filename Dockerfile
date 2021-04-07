FROM nfcore/base:1.13.3
LABEL authors="Yanhai Gong" \
      description="Docker image containing all requirements for gongyh/nf-core-scgs pipeline"

# Install procps so that Nextflow can poll CPU usage, set locale en_US.UTF-8 used by Picard
RUN apt-get update && apt-get install --no-install-recommends -y procps libpng16-16 cmake gcc g++ hmmer2 locales && \
  apt-get autoremove --purge && apt-get clean -y && apt-get autoremove && rm -rf /var/lib/apt/lists/* && \
  sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen && locale-gen

# Install conda environments
RUN conda install -y mamba nomkl -c conda-forge && conda clean -y -a && rm -rf /opt/conda/pkgs/*
COPY environment.yml /

RUN mamba env create --quiet -f /environment.yml && conda clean -a
RUN echo 'conda activate nf-core-gongyh-scgs-1.1.3' >> ~/.bashrc
ENV PATH /opt/conda/envs/nf-core-gongyh-scgs-1.1.3/bin:$PATH

#RUN mamba env create -n scgs_py36 -f /environment.yml && \
#    conda clean -y -a && rm -rf /opt/conda/pkgs/*
#RUN echo 'conda activate scgs_py36' >> ~/.bashrc
#ENV PATH /opt/conda/envs/scgs_py36/bin:$PATH

# Install Bioconductor packages
RUN R -e "install.packages(c('BiocManager','stringi'), repos='https://cloud.r-project.org'); BiocManager::install('GenomeInfoDbData'); BiocManager::install('AneuFinder')"

# Download silva and busco for Quast 5.x, update taxa db for krona
RUN quast-download-silva && ktUpdateTaxonomy.sh #&& quast-download-busco

# Download GATK3 and setting
#RUN wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-0-ge9d806836" -O GenomeAnalysisTK-3.8.tar.bz2 && \
#    tar xjvf GenomeAnalysisTK-3.8.tar.bz2 && gatk3-register GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar && rm -rf ./GenomeAnalysisTK-3.8*
COPY GenomeAnalysisTK.jar /
RUN gatk3-register /GenomeAnalysisTK.jar && rm -rf /GenomeAnalysisTK.jar

# Install conda environment py27
COPY py27_env.yml /
RUN mamba env create -n scgs_py27 -f /py27_env.yml nomkl && conda clean -y -a && rm -rf /opt/conda/pkgs/*

# Install ACDC
RUN git clone https://github.com/mlux86/acdc.git /tmp/acdc && cd /tmp/acdc && mkdir build && cd build && cmake .. -DBOOST_ROOT=/opt/conda/envs/nf-core-gongyh-scgs-1.1.3/ && \
    make -j $(nproc) && make install && rm -rf /tmp/acdc

# Install rnammer
RUN mkdir /tmp/rnammer
ADD rnammer-1.2.src.tar.Z /tmp/rnammer
ADD rnammer.patch /tmp/rnammer
RUN cd /tmp/rnammer && tar xf rnammer-1.2.src.tar.Z && rm rnammer-1.2.src.tar.Z && patch < rnammer.patch && mkdir /usr/local/share/rnammer && \
    cp -r * /usr/local/share/rnammer && ln -s /usr/local/share/rnammer/rnammer /usr/local/bin/rnammer && rm -rf /tmp/rnammer

# Build database for blobtools
RUN [ "/bin/bash", "-c", "source activate scgs_py27 && blobtools-build_nodesdb && source deactivate" ]

# Install Funannotate
#RUN [ "/bin/bash", "-c", "source activate scgs_py27 && mamba install -y -c bioconda -c conda-forge trnascan-se=1.3.1 funannotate=1.7.2 && conda clean -y -a && rm -rf /opt/conda/pkgs/* && source deactivate" ]

# Install R packages
RUN R -e "install.packages(c('magicaxis','ape','gridExtra','genoPlotR','ggpubr','rstatix'), repos='https://cloud.r-project.org')"

# Install packages
RUN apt update && apt install -y autoconf automake unzip zlib1g-dev libncurses5-dev libncursesw5-dev fastx-toolkit \
    augustus augustus-data augustus-doc && apt-get autoremove --purge && apt-get clean && apt-get autoremove

# Install ANIcalculator
ADD ANIcalculator_v1.tgz /usr/local/bin/
RUN cd /usr/local/bin && cp ANIcalculator_v1/ANIcalculator . && cp ANIcalculator_v1/nsimscan . && rm -rf ANIcalculator_v1

# Install packages from git source
RUN cd /opt && git clone https://github.com/gongyh/ezTree.git && git clone https://bitbucket.org/genomicepidemiology/virulencefinder.git && \
    git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git && git clone https://bitbucket.org/genomicepidemiology/pointfinder.git

# Keep a copy of current pipeline to container
COPY . /opt/nf-core-scgs/

# store conda environments
RUN conda env export --name nf-core-gongyh-scgs-1.1.3 > nf-core-gongyh-scgs-1.1.3.yml
RUN conda env export --name scgs-py27 > scgs-py27.yml

# Copy conda (de)activate script to global PATH
RUN cp /opt/conda/bin/activate /usr/local/bin/ && cp /opt/conda/bin/deactivate /usr/local/bin/

# Add default container command
CMD ["/bin/bash"]

