FROM ubuntu:22.04

LABEL maintainer="joao.cardoso@example.com"
LABEL description="Cell Ranger 10.0.0 for 10x Genomics Perturb-Seq analysis (industry standard)"
LABEL version="1.0"

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    tar \
    unzip \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Download and install Cell Ranger 10.0.0
RUN wget -O /tmp/cellranger-10.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1768369176&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=m5FIq9csLJgXUnnRR1jCOez-nbmt83shOTOSzGTeAs1o8IMgHw8w1xjg3tP6P~wTtel6tVEMGz75q554YCqR-Qk5BoMfwxpm7YFBvRqYU-d79Be1CdVhUbILHqQq8Ip32wS0M6weD2e06YiuiR2GgQDQc3Ss8nvj7cPEKv92pWCFZB9lm5jFYFFYMQ6WBYmqnEaTbfm-HtrE4q6c4GFLqTQtKlRZQqM9W8aLHDE6SCqgepPJZV9MsJCaTdthJHtHLPawUAQJ-eKNHyePT-ZWzCZJWpRRxOfog71E5xlJrrrWERYdPMnx~VtDE-eeT1BEpVPz1MpqoL9S5oLBQigkAA__" && \
    tar -xzf /tmp/cellranger-10.0.0.tar.gz -C /opt/ && \
    rm /tmp/cellranger-10.0.0.tar.gz

# Add Cell Ranger to PATH
ENV PATH="/opt/cellranger-10.0.0:${PATH}"

# Note: Verification is skipped during build to avoid Rosetta issues on ARM64
# Cell Ranger will be verified when the container runs

WORKDIR /data

CMD ["bash"]
