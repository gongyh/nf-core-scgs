FROM nfcore/base
LABEL authors="Yanhai Gong" \
      description="Docker image containing all requirements for nf-core/scgs pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-scgs-1.0dev/bin:$PATH
