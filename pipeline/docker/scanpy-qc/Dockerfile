FROM python:3.10-slim

LABEL maintainer="joao.cardoso@example.com"
LABEL description="Scanpy for single-cell RNA-seq QC and filtering"
LABEL version="1.0"

# Install system dependencies for scientific computing
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    libhdf5-dev \
    pkg-config \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for single-cell analysis
RUN pip install --no-cache-dir \
    scanpy==1.9.6 \
    pandas==2.1.4 \
    numpy==1.26.2 \
    matplotlib==3.8.2 \
    scipy==1.11.4 \
    scikit-learn==1.3.2 \
    scikit-misc==0.3.1 \
    anndata==0.10.3 \
    leidenalg \
    python-igraph \
    tables==3.9.2 \
    h5py==3.10.0

# Create directory for scripts
RUN mkdir -p /scripts
WORKDIR /scripts

# Copy QC filter script
COPY scripts/qc_filter.py /scripts/

# Set working directory for data processing
WORKDIR /data

# Test that packages are installed
RUN python -c "import scanpy; import anndata; print(f'Scanpy {scanpy.__version__} installed successfully')"

CMD ["python", "--version"]
