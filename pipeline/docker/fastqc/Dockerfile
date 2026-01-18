FROM ubuntu:22.04

LABEL maintainer="joao.cardoso@example.com"
LABEL description="FastQC for quality control of FASTQ sequencing data"
LABEL version="1.0"

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    openjdk-11-jre \
    perl \
    fontconfig \
    libfindbin-libs-perl \
    && rm -rf /var/lib/apt/lists/*

# Install FastQC
WORKDIR /opt
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod +x FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.12.1.zip

# Set working directory for data processing
WORKDIR /data

# Test that FastQC is installed correctly
RUN fastqc --version

# Default command shows help
CMD ["fastqc", "--help"]
