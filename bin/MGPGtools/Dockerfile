FROM mambaorg/micromamba:latest

LABEL maintainer "Yanhai Gong <gongyh@qibebt.ac.cn>"

COPY /environment.yml /

RUN micromamba create -f /environment.yml && micromamba clean -a

COPY / .

RUN eval "$(micromamba shell hook --shell=bash)" && micromamba activate mgpg-tools && python setup.py install

ENV ENV_NAME=mgpg-tools
