FROM python:3.10-slim

LABEL maintainer="joao.cardoso@example.com"
LABEL description="Single-cell analysis: clustering, dimensionality reduction, and visualization"
LABEL version="1.0"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    libhdf5-dev \
    pkg-config \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install comprehensive Python packages for analysis
RUN pip install --no-cache-dir \
    scanpy==1.9.6 \
    pandas==2.1.4 \
    numpy==1.26.2 \
    matplotlib==3.8.2 \
    scipy==1.11.4 \
    scikit-learn==1.3.2 \
    anndata==0.10.3 \
    leidenalg \
    python-igraph \
    umap-learn==0.5.5 \
    plotly==5.18.0 \
    tables==3.9.2 \
    h5py==3.10.0 \
    adjustText==0.8

# Create directory for scripts
RUN mkdir -p /scripts
WORKDIR /scripts

# Scripts will be copied during build
COPY scripts/clustering.py /scripts/

# Set working directory
WORKDIR /data

# Test installations
RUN python -c "import scanpy; import umap; import leidenalg; print('All packages installed successfully')"

CMD ["python", "--version"]
