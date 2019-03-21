FROM continuumio/miniconda3:4.5.4
LABEL authors="Yanhai Gong" \
      description="Docker image containing all requirements for gongyh/nf-core-scgs pipeline"

COPY environment.yml /
RUN conda env update -n base -f /environment.yml && conda clean -a

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y 

# Install graphviz so that Nextflow can draw DAG
RUN apt-get update && apt-get install -y graphviz && apt-get clean

# Install Circleator
RUN apt-get update && apt-get install -y libbatik-java vcftools && apt-get clean
RUN cpanm JSON Log::Log4perl SVG Text::CSV Test::Most Bio::FeatureIO::gff
RUN cd /opt && wget https://github.com/jonathancrabtree/Circleator/archive/1.0.1.tar.gz -O Circleator-1.0.1.tar.gz && \
    tar xzvf Circleator-1.0.1.tar.gz && cd Circleator-1.0.1 && perl Build.PL && ./Build && ./Build install
RUN rm -rf /opt/*.tar.gz

# clean up
RUN apt-get autoremove --purge && apt-get clean && apt-get autoremove && conda clean -a

