# Força arquitetura compatível com bcl2fastq2 (x86_64)
FROM mambaorg/micromamba:latest

LABEL maintainer="joao.cardoso@example.com"
LABEL description="Illumina bcl2fastq2 for BCL to FASTQ conversion (primary analysis)"
LABEL version="1.0"

USER root

# Canais + prioridade estrita (boa prática)
RUN micromamba config append channels conda-forge && \
    micromamba config append channels bioconda && \
    micromamba config append channels defaults && \
    micromamba config append channels bih-cubi && \
    micromamba config set channel_priority strict

# Instalação (podes pinhar versão)
RUN micromamba install -y -n base bcl2fastq2=2.20.0.422 procps-ng && \
    micromamba clean -a -y

# Verificação
RUN micromamba run -n base bcl2fastq --version

WORKDIR /data

# Recomendação do micromamba-docker: garantir que os comandos correm no env
CMD ["bash"]