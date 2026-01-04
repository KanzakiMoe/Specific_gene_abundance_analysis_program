FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
SHELL ["/bin/bash", "-c"]

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    python3 \
    python3-pip \
    samtools \
    bowtie2 \
    hmmer \
    prodigal \
    fastp \
    zlib1g-dev \
    megahit \
    dos2unix \
&& rm -rf /var/lib/apt/lists/*

# Miniconda (for install CoverM due to apt-get version is unable to use)
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/opt/conda/bin:$PATH

# CoverM (due to apt-get version is unable to use)
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda install -y coverm && \
    conda clean -afy


# work path definition
WORKDIR /work

# pipeline script
COPY run_pipeline.sh /work/run_pipeline.sh
COPY plot_anammox_abundance.R /work/plot_anammox_abundance.R
RUN chmod +x /work/run_pipeline.sh

ENTRYPOINT ["/bin/bash", "-c", "\
    find /work -type f -name '*.sh' -exec dos2unix {} \\; || true; \
    bash /work/run_pipeline.sh \
"]

# R for plotting
RUN apt-get install -y r-base \
    && R -e "install.packages(c('stringr','ggplot2','dplyr','readr','tidyr'), repos='https://cloud.r-project.org')"

