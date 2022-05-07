FROM nfcore/base:1.13.3
LABEL authors="Yanhai Gong" \
      description="Docker image containing all requirements for gongyh/nf-core-scgs pipeline"

# Install procps so that Nextflow can poll CPU usage, set locale en_US.UTF-8 used by Picard
RUN apt-get update && apt-get install --no-install-recommends -y procps libpng16-16 \
    cmake gcc g++ hmmer2 locales autoconf automake unzip zlib1g-dev libncurses5-dev \
    libncursesw5-dev fastx-toolkit augustus augustus-data augustus-doc libidn11 libgl1 && \
  apt-get autoremove --purge && apt-get clean -y && apt-get autoremove && rm -rf /var/lib/apt/lists/* && \
  sed -i 's/^# *\(en_US.UTF-8\)/\1/' /etc/locale.gen && locale-gen

# Install conda environments
RUN conda install -y mamba nomkl eggnog-mapper=2.1 -c conda-forge -c bioconda && conda clean -y -a && rm -rf /opt/conda/pkgs/*
COPY environment.yml /
RUN mamba env create --quiet -f /environment.yml && conda clean -a
RUN echo 'conda activate nf-core-gongyh-scgs' >> ~/.bashrc
ENV PATH /opt/conda/envs/nf-core-gongyh-scgs/bin:$PATH

# Download silva and busco for Quast 5.x, update taxa db for krona
RUN quast-download-silva && ktUpdateTaxonomy.sh && \
    sed -i 's#/bin/blobtools#/lib/python3.6/site-packages#g' /opt/conda/envs/nf-core-gongyh-scgs/bin/blobtools-build_nodesdb && \
    sed -i 's#$blobtools/lib/blobtools.py#blobtools#g' /opt/conda/envs/nf-core-gongyh-scgs/bin/blobtools-build_nodesdb && \
    blobtools-build_nodesdb

# Install ACDC
RUN git config --global http.sslVerify false && git clone https://github.com/mlux86/acdc.git /tmp/acdc && \
    cd /tmp/acdc && mkdir build && cd build && \
    cmake -DBoost_NO_BOOST_CMAKE=TRUE -DBOOST_ROOT:PATHNAME=/opt/conda/envs/nf-core-gongyh-scgs/ -DBoost_NO_SYSTEM_PATHS=ON .. && \
    make -j $(nproc) && make install && rm -rf /tmp/acdc

# Install rnammer
RUN mkdir /tmp/rnammer
ADD rnammer-1.2.src.tar.Z /tmp/rnammer
ADD rnammer.patch /tmp/rnammer
RUN cd /tmp/rnammer && tar xf rnammer-1.2.src.tar.Z && rm rnammer-1.2.src.tar.Z && patch < rnammer.patch && mkdir /usr/local/share/rnammer && \
    cp -r * /usr/local/share/rnammer && ln -s /usr/local/share/rnammer/rnammer /usr/local/bin/rnammer && rm -rf /tmp/rnammer

# Install ANIcalculator
ADD ANIcalculator_v1.tgz /usr/local/bin/
RUN cd /usr/local/bin && cp ANIcalculator_v1/ANIcalculator . && cp ANIcalculator_v1/nsimscan . && rm -rf ANIcalculator_v1

# Install packages from git source
RUN cd /opt && git config --global http.sslVerify false && \
    git clone https://github.com/gongyh/ezTree.git && \
    git clone -b 2.0.4 https://bitbucket.org/genomicepidemiology/virulencefinder.git && \
    git clone -b 3.2.1 https://git@bitbucket.org/genomicepidemiology/resfinder.git && \
    git clone -b 3.1.1 https://bitbucket.org/genomicepidemiology/pointfinder.git

# Install Monovar
RUN cd /opt && git config --global http.sslVerify false && \
    git clone https://github.com/gongyh/MonoVar.git && \
    cd MonoVar && python setup.py install

# Keep a copy of current pipeline to container
COPY . /opt/nf-core-scgs/

# store conda environments
RUN conda env export --name nf-core-gongyh-scgs > nf-core-gongyh-scgs.yml

# Copy conda (de)activate script to global PATH
RUN cp /opt/conda/bin/activate /usr/local/bin/ && cp /opt/conda/bin/deactivate /usr/local/bin/

# Add default container command
CMD ["/bin/bash"]

