# syntax=docker/dockerfile:1
ARG R_VERSION=4.5.2
FROM rocker/rstudio:${R_VERSION}

ARG NCPUS=4
ENV DEBIAN_FRONTEND=noninteractive \
    BIOCONDUCTOR_VERSION=3.22 \
    CRAN_REPO=https://p3m.dev/cran/__linux__/noble/2026-03-10 \
    MAKEFLAGS=-j2 \
    RGL_USE_NULL=TRUE \
    TZ=Europe/Rome \
    RETICULATE_PYTHON=/opt/catcher-python/bin/python3
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash coreutils sed grep gawk gzip bzip2 pigz unzip zip tar \
    git curl wget ca-certificates python3 python3-dev python3-venv python3-pip \
    build-essential gfortran gobjc++ cmake pkg-config \
    zlib1g-dev libbz2-dev liblzma-dev libpcre2-dev libreadline-dev \
    libcurl4-openssl-dev libxml2-dev libssl-dev libgit2-dev \
    libudunits2-dev libgdal-dev libgeos-dev libproj-dev \
    libharfbuzz-dev libfribidi-dev libcairo2-dev libfontconfig1-dev \
    libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev \
    libglpk-dev libhdf5-dev libfftw3-dev libgmp3-dev libmpfr-dev \
    libtbb-dev libblas-dev liblapack-dev libnlopt-dev libgsl-dev \
    libicu-dev libsqlite3-dev libzstd-dev libx11-dev libxt-dev \
    libgl1-mesa-dev libglu1-mesa-dev fonts-liberation \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m venv /opt/catcher-python
ENV PATH="/opt/catcher-python/bin:${PATH}"

COPY packages-*.txt sources.json /opt/catcher-build/
COPY scripts/install-r.R scripts/verify-image.R /opt/catcher-build/
RUN NCPUS=${NCPUS} Rscript /opt/catcher-build/install-r.R cran
# Required by the CRAN magick binary, imported by SpatialExperiment/Banksy.
RUN apt-get update && apt-get install -y --no-install-recommends libmagick++-dev \
    && rm -rf /var/lib/apt/lists/*
RUN NCPUS=${NCPUS} Rscript /opt/catcher-build/install-r.R bioc
RUN --mount=type=bind,source=vendor,target=/opt/catcher-vendor,readonly \
    NCPUS=${NCPUS} Rscript /opt/catcher-build/install-r.R github \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

# tidyseurat e monocle3: installati per ultimi per non rompere le dipendenze
RUN Rscript -e 'install.packages("tidyseurat")' \
    && Rscript -e 'devtools::install_github("cole-trapnell-lab/monocle3")' \
    && rm -rf /tmp/downloaded_packages /tmp/Rtmp*

RUN Rscript /opt/catcher-build/verify-image.R \
    && python3 -m compileall -q /usr/local/lib/R/site-library/catcheR/python \
    && bash --version > /opt/catcher-build/bash-version.txt \
    && python3 --version > /opt/catcher-build/python-version.txt \
    && ! command -v cellranger \
    && ! command -v bcl2fastq

RUN mkdir -p /home/rstudio/project /opt/catcher/templates \
    && cp -a /usr/local/lib/R/site-library/catcheR/rmd/. /opt/catcher/templates/ \
    && mkdir -p /data \
    && ln -s /home/rstudio/project/scratch /data/scratch \
    && for f in /usr/local/lib/R/site-library/catcheR/scripts/*; do ln -s "$f" /home/"$(basename "$f")"; done \
    && chmod +x /usr/local/lib/R/site-library/catcheR/scripts/*.sh \
    && chown -R rstudio:rstudio /home/rstudio/project
COPY catcheR_V4/ /opt/catcheR/source/
RUN ln -s /opt/catcheR/source /home/rstudio/catcheR_V4

LABEL org.opencontainers.image.title="catcheR V4 RStudio" \
      org.opencontainers.image.description="RStudio, full analysis R library collection, Python and Bash" \
      org.opencontainers.image.source="https://github.com/ebalmas/catcheR" \
      org.opencontainers.image.revision="724ef8900a6edb31136a7b0ccf4d3d5e3ba18dbe"
WORKDIR /home/rstudio/project
EXPOSE 8787
CMD ["/init"]
