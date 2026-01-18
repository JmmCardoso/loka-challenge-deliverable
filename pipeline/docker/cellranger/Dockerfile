FROM ubuntu:20.04

LABEL maintainer="GeneXOmics Pipeline"
LABEL description="Cell Ranger 10.0.0 for 10x Genomics single-cell analysis"

# Prevent interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    tar \
    gzip \
    debianutils \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Set Cell Ranger version
ENV CELLRANGER_VER=10.0.0

# Download Cell Ranger directly from 10x Genomics
# Note: If download fails, visit https://www.10xgenomics.com/support/software/cell-ranger/downloads
# to get a fresh download link
RUN cd /opt && \
    curl -o cellranger-${CELLRANGER_VER}.tar.gz \
    "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1768797903&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=D-J~~83kXIODkZ25nFjahASWAdYi-nV7P6WkSPp-ceymBNI7h3GKVWjPrJsuRbtHposa4ru9EGShX-h9xyM8seeVZSYocYne~k45zAw0AqaGOA-biNItc7BYjaqaPpMaMcT8Icf7zVa5GQaVO6qrLRM7J28w6DvQaVGULLqGthWADbJ76fI~4g1Y6Z2bDhdXyZ9VTJ8fr6o91w4150zJpoZE1oJ69jqP5cH1VzPFiiBLL8rRzXfX23AFIEedBMGpefeflZWKEfRkTbJyG37sZs2Wdq1pIX7xZlZUbk5-dsQBvRBTbBPE~dpwOhUv6mlD0Fb3DxfwX-3doWeutwWSBw__" && \
    tar -xzf cellranger-${CELLRANGER_VER}.tar.gz && \
    rm cellranger-${CELLRANGER_VER}.tar.gz

# Add Cell Ranger to PATH
ENV PATH="/opt/cellranger-${CELLRANGER_VER}:${PATH}"

# Create symbolic link for easier access
RUN ln -s /opt/cellranger-${CELLRANGER_VER}/cellranger /usr/local/bin/cellranger

# Note: Verification skipped due to ARM64/x86-64 compatibility
# Cell Ranger will run when executed by Nextflow with proper platform settings

# Set working directory
WORKDIR /data

# Default command
CMD ["cellranger", "--help"]
