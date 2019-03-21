FROM nfcore/base
LABEL authors="Yanhai Gong" \
      description="Docker image containing all requirements for gongyh/nf-core-scgs pipeline"

RUN apt-get update && apt-get install -y r-base && apt-get clean

RUN apt-get update && apt-get install -y graphviz && apt-get clean

RUN apt-get update && apt-get install -y perl bioperl cpanminus libbatik-java vcftools && apt-get clean
RUN cpanm JSON Log::Log4perl SVG Text::CSV Test::Most Bio::FeatureIO::gff
RUN cd /opt && wget https://github.com/jonathancrabtree/Circleator/archive/1.0.1.tar.gz -O Circleator-1.0.1.tar.gz && \
    tar xzvf Circleator-1.0.1.tar.gz && cd Circleator-1.0.1 && perl Build.PL && ./Build && ./Build install
RUN rm -rf /opt/*.tar.gz

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

RUN apt-get autoremove --purge && apt-get clean && apt-get autoremove && conda clean -a

ENTRYPOINT conda activate nf-core-scgs-1.0dev
ENV PATH /opt/conda/envs/nf-core-scgs-1.0dev/bin:$PATH
