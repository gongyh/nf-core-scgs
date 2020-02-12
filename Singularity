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
    apt-get update && apt-get install -y procps libpng16-16 cmake gcc g++ hmmer2 && apt-get clean -y
    apt-get install -y locales && sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen && locale-gen
    conda install python=3.6 conda=4.6.12 && conda env update -n base -f /environment.yml && conda clean -a
    R -e "install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('GenomeInfoDbData'); BiocManager::install('AneuFinder')"
    quast-download-silva && ktUpdateTaxonomy.sh
    #wget "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-0-ge9d806836" -O GenomeAnalysisTK-3.8.tar.bz2 && \
    #  tar xjvf GenomeAnalysisTK-3.8.tar.bz2 && gatk3-register GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar && rm -rf ./GenomeAnalysisTK-3.8*
    gatk3-register /opt/nf-core-scgs/GenomeAnalysisTK.jar
    conda env create -n py27 -f /py27_env.yml && conda clean -a
    /bin/bash -c "source activate py27 && blobtools-build_nodesdb && source deactivate"
    /bin/bash -c "source activate py27 && conda install -c bioconda -c conda-forge funannotate=1.7.2 && source deactivate"
    cd /opt && git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git
    cd /opt && git clone https://bitbucket.org/genomicepidemiology/pointfinder.git
    git clone https://github.com/mlux86/acdc.git /tmp/acdc && cd /tmp/acdc && mkdir build && cd build && cmake .. -DBOOST_ROOT=/opt/conda/ && \
      make -j $(nproc) && make install && rm -rf /tmp/acdc && mkdir /acdc
    mkdir /tmp/rnammer && cd /tmp/rnammer && cp /opt/nf-core-scgs/rnammer-1.2.src.tar.Z . && cp /opt/nf-core-scgs/rnammer.patch .
    tar xf rnammer-1.2.src.tar.Z && rm rnammer-1.2.src.tar.Z && patch < rnammer.patch && mkdir /usr/local/share/rnammer && \
      cp -r * /usr/local/share/rnammer && ln -s /usr/local/share/rnammer/rnammer /usr/local/bin/rnammer && rm -rf /tmp/rnammer
    cd /opt && wget https://github.com/takaram/kofam_scan/archive/v1.1.0.tar.gz -O kofamscan-1.1.0.tar.gz && \
      tar --no-same-owner -xzvf kofamscan-1.1.0.tar.gz && rm -rf kofamscan-1.1.0.tar.gz
    apt-get autoremove --purge && apt-get clean && apt-get autoremove
    conda clean -y -a && rm -rf /opt/conda/pkgs/*
    chmod -R o+rx /opt/nf-core-scgs

%runscript
    exec /opt/conda/bin/nextflow run /opt/nf-core-scgs/main.nf --help
