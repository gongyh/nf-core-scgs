Bootstrap: docker
FROM: continuumio/miniconda3:4.5.4
IncludeCmd: yes

%help
    Singularity container for gongyh/scgs pipeline.

%files
    environment.yml /
    py27_env.yml /
    . /opt/nf-core-scgs/

%labels
    Author gongyh@qibebt.ac.cn

%post
    export PATH="/opt/conda/bin:$PATH"
    apt-get update && apt-get install --no-install-recommends -y procps libpng16-16 cmake gcc g++ hmmer2 locales && \
      apt-get autoremove --purge && apt-get clean -y && apt-get autoremove && rm -rf /var/lib/apt/lists/* && \
      sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen && locale-gen
    cd /opt && git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
    cd /opt && git clone https://bitbucket.org/genomicepidemiology/pointfinder.git
    conda install python=3.6 conda=4.7.12 nomkl && conda install -y mamba -c conda-forge && conda clean -y -a && rm -rf /opt/conda/pkgs/*
    mamba env create -n scgs_py36 -f /environment.yml && conda clean -y -a && rm -rf /opt/conda/pkgs/*
    echo 'conda activate scgs_py36' >> /root/.bashrc
    export PATH="/opt/conda/envs/scgs_py36/bin:$PATH"
    R -e "install.packages(c('BiocManager','stringi'), repos='https://cloud.r-project.org'); BiocManager::install('GenomeInfoDbData'); BiocManager::install('AneuFinder')"
    quast-download-silva && ktUpdateTaxonomy.sh
    #wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-0-ge9d806836" -O GenomeAnalysisTK-3.8.tar.bz2 && \
    #  tar xjvf GenomeAnalysisTK-3.8.tar.bz2 && gatk3-register GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar && rm -rf ./GenomeAnalysisTK-3.8*
    #gatk3-register /opt/nf-core-scgs/GenomeAnalysisTK.jar
    mamba env create -n scgs_py27 -f /py27_env.yml nomkl && conda clean -y -a && rm -rf /opt/conda/pkgs/*
    git clone https://github.com/mlux86/acdc.git /tmp/acdc && cd /tmp/acdc && mkdir build && cd build && cmake .. -DBOOST_ROOT=/opt/conda/envs/scgs_py36/ && \
      make -j $(nproc) && make install && rm -rf /tmp/acdc && mkdir /acdc
    mkdir /tmp/rnammer && cd /tmp/rnammer && cp /opt/nf-core-scgs/rnammer-1.2.src.tar.Z . && cp /opt/nf-core-scgs/rnammer.patch . && \
      tar xf rnammer-1.2.src.tar.Z && rm rnammer-1.2.src.tar.Z && patch < rnammer.patch && mkdir /usr/local/share/rnammer && \
      cp -r * /usr/local/share/rnammer && ln -s /usr/local/share/rnammer/rnammer /usr/local/bin/rnammer && rm -rf /tmp/rnammer
    /bin/bash -c "source activate scgs_py27 && blobtools-build_nodesdb && source deactivate"
    /bin/bash -c "source activate scgs_py27 && mamba install -y -c bioconda -c conda-forge trnascan-se=1.3.1 funannotate=1.7.2 && conda clean -y -a && rm -rf /opt/conda/pkgs/* && source deactivate"
    cd /opt && wget https://github.com/takaram/kofam_scan/archive/v1.1.0.tar.gz -O kofamscan-1.1.0.tar.gz && \
      tar --no-same-owner -xzvf kofamscan-1.1.0.tar.gz && rm -rf kofamscan-1.1.0.tar.gz
    R -e "install.packages(c('magicaxis','ape','gridExtra','hyperSpec','permute'), repos='https://cloud.r-project.org')"
    cd /opt && git clone --recursive https://github.com/mcveanlab/mccortex && cd mccortex && \
      apt update && apt install -y autoconf automake zlib1g-dev libncurses5-dev libncursesw5-dev fastx-toolkit && \
      make all && apt-get autoremove --purge && apt-get clean && apt-get autoremove
    cd /usr/local/bin && wget http://opengene.org/fastp/fastp && chmod a+x ./fastp
    conda clean -y -a && rm -rf /opt/conda/pkgs/*
    chmod -R o+rx /opt/nf-core-scgs

%runscript
    exec /bin/bash
