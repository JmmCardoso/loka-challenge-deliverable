#!/bin/bash

#===============================================================================
# Setup Test Data for Perturb-Seq Pipeline
#===============================================================================
# Downloads all required test datasets and reference genome for pipeline testing:
# 1. 10x Genomics Perturb-Seq test data (CRISPR + GEX)
# 2. MiSeq ATAC BCL test data (for BCL conversion testing)
# 3. Cell Ranger reference genome (GRCh38-2024-A)
#
# Author: João Cardoso
# Date: January 2026
#===============================================================================

set -e  # Exit on error

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="${PIPELINE_DIR}/data"

echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}  Perturb-Seq Pipeline - Test Data Setup${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""
echo "This script will download:"
echo "  1. 10x Genomics Perturb-Seq test data (~5 GB)"
echo "  2. MiSeq ATAC BCL test data (~500 MB)"
echo "  3. Cell Ranger reference genome (~15 GB)"
echo ""
echo "Total download: ~20 GB"
echo "Estimated time: 30-60 minutes (depending on connection)"
echo ""

# Ask for confirmation
printf "${YELLOW}Continue with download? (y/n) ${NC}"
read -r REPLY
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Setup cancelled."
    exit 0
fi

# Create data directories
echo -e "${GREEN}Creating directory structure in ${DATA_DIR}...${NC}"
mkdir -p "${DATA_DIR}/10x_test_data"
mkdir -p "${DATA_DIR}/miseq-atac"
mkdir -p "${DATA_DIR}/references"

#===============================================================================
# 1. Download 10x Genomics Perturb-Seq Test Data
#===============================================================================

echo ""
echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}[1/3] Downloading 10x Genomics Perturb-Seq Test Data${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""

if [ -d "${DATA_DIR}/10x_test_data/1k_CRISPR_5p_gemx_fastqs" ]; then
    echo -e "${YELLOW}10x test data already exists. Skipping download.${NC}"
else
    echo "Dataset: 1k CRISPR 5' Gene Expression Multiplex"
    echo "Source: 10x Genomics Public Datasets"
    echo "Size: ~5 GB"
    echo ""
    
    cd "${DATA_DIR}/10x_test_data"
    
    # Download from 10x Genomics
    echo "Downloading FASTQ files..."
    wget -O 1k_CRISPR_5p_gemx_fastqs.tar \
        "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/1k_CRISPR_5p_gemx_Multiplex/1k_CRISPR_5p_gemx_Multiplex_fastqs.tar" || {
        echo -e "${RED}Error: Failed to download 10x test data${NC}"
        exit 1
    }
    
    echo "Extracting files..."
    tar -xf 1k_CRISPR_5p_gemx_fastqs.tar
    rm 1k_CRISPR_5p_gemx_fastqs.tar
    
    # Download feature reference and config
    echo "Downloading feature reference..."
    wget -O 1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv \
        "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/1k_CRISPR_5p_gemx_Multiplex/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv"
    
    wget -O 1k_CRISPR_5p_gemx_Multiplex_config.csv \
        "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/1k_CRISPR_5p_gemx_Multiplex/1k_CRISPR_5p_gemx_Multiplex_config.csv"
    
    cd "${PIPELINE_DIR}"
    
    echo -e "${GREEN}✓ 10x test data download complete${NC}"
fi

#===============================================================================
# 2. Download MiSeq ATAC BCL Test Data
#===============================================================================

echo ""
echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}[2/3] Downloading MiSeq ATAC BCL Test Data${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""

if [ -f "${DATA_DIR}/miseq-atac/RunInfo.xml" ]; then
    echo -e "${YELLOW}MiSeq BCL data already exists. Skipping download.${NC}"
else
    echo "Dataset: MiSeq ATAC-seq BCL files"
    echo "Source: Illumina BaseSpace (sample data)"
    echo "Size: ~500 MB"
    echo ""
    
    cd "${DATA_DIR}"
    
    # Download Cell Ranger ARC tiny BCL for ATAC
    echo "Downloading BCL files..."
    wget -O miseq-atac-bcl.tar.gz \
        "https://cf.10xgenomics.com/supp/cell-arc/cellranger-arc-tiny-bcl-atac-2.0.0.tar.gz" || {
        echo -e "${RED}Error: Failed to download MiSeq BCL data${NC}"
        exit 1
    }
    
    echo "Extracting files..."
    tar -xzf miseq-atac-bcl.tar.gz
    rm miseq-atac-bcl.tar.gz
    
    # Create compatible sample sheet if needed
    if [ ! -f "${DATA_DIR}/miseq-atac/SampleSheet.csv" ]; then
        cp "${DATA_DIR}/miseq-atac/IlluminaSampleSheet.csv" "${DATA_DIR}/miseq-atac/SampleSheet.csv"
    fi
    
    cd "${PIPELINE_DIR}"
    
    echo -e "${GREEN}✓ MiSeq BCL data download complete${NC}"
fi

#===============================================================================
# 3. Download Cell Ranger Reference Genome
#===============================================================================

echo ""
echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}[3/3] Downloading Cell Ranger Reference Genome${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""

if [ -d "${DATA_DIR}/references/refdata-gex-GRCh38-2024-A" ]; then
    echo -e "${YELLOW}Reference genome already exists. Skipping download.${NC}"
else
    echo "Reference: Human GRCh38 (2024-A)"
    echo "Source: 10x Genomics"
    echo "Size: ~15 GB"
    echo ""
    echo -e "${YELLOW}This is the largest file and will take time...${NC}"
    echo ""
    
    cd "${DATA_DIR}/references"
    
    # Download reference genome
    echo "Downloading Human reference genome..."
    wget -O refdata-gex-GRCh38-2024-A.tar.gz \
        "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz" || {
        echo -e "${RED}Error: Failed to download reference genome${NC}"
        echo "You can manually download from: https://www.10xgenomics.com/support/software/cell-ranger/downloads"
        exit 1
    }
    
    echo "Extracting reference genome (this may take 5-10 minutes)..."
    tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
    rm refdata-gex-GRCh38-2024-A.tar.gz
    
    cd "${PIPELINE_DIR}"
    
    echo -e "${GREEN}✓ Reference genome download complete${NC}"
fi

#===============================================================================
# Verification and Summary
#===============================================================================

echo ""
echo -e "${BLUE}================================================================${NC}"
echo -e "${BLUE}  Setup Complete - Data Summary${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""

# Check 10x data
if [ -d "${DATA_DIR}/10x_test_data/1k_CRISPR_5p_gemx_fastqs" ]; then
    echo -e "${GREEN}✓${NC} 10x Genomics Perturb-Seq data"
    echo "   Location: data/10x_test_data/"
    CRISPR_COUNT=$(find "${DATA_DIR}/10x_test_data/1k_CRISPR_5p_gemx_fastqs/crispr" -name "*.fastq.gz" 2>/dev/null | wc -l)
    GEX_COUNT=$(find "${DATA_DIR}/10x_test_data/1k_CRISPR_5p_gemx_fastqs/gex" -name "*.fastq.gz" 2>/dev/null | wc -l)
    echo "   Files: ${CRISPR_COUNT} CRISPR FASTQs, ${GEX_COUNT} GEX FASTQs"
else
    echo -e "${RED}✗${NC} 10x test data missing"
fi

echo ""

# Check BCL data
if [ -f "${DATA_DIR}/miseq-atac/RunInfo.xml" ]; then
    echo -e "${GREEN}✓${NC} MiSeq ATAC BCL data"
    echo "   Location: data/miseq-atac/"
    BCL_SIZE=$(du -sh "${DATA_DIR}/miseq-atac" 2>/dev/null | cut -f1)
    echo "   Size: ${BCL_SIZE}"
else
    echo -e "${RED}✗${NC} MiSeq BCL data missing"
fi

echo ""

# Check reference
if [ -d "${DATA_DIR}/references/refdata-gex-GRCh38-2024-A" ]; then
    echo -e "${GREEN}✓${NC} Cell Ranger reference genome"
    echo "   Location: data/references/refdata-gex-GRCh38-2024-A/"
    REF_SIZE=$(du -sh "${DATA_DIR}/references/refdata-gex-GRCh38-2024-A" 2>/dev/null | cut -f1)
    echo "   Size: ${REF_SIZE}"
else
    echo -e "${RED}✗${NC} Reference genome missing"
fi

echo ""
echo -e "${BLUE}================================================================${NC}"
echo -e "${GREEN}All test data successfully downloaded!${NC}"
echo -e "${BLUE}================================================================${NC}"
echo ""
echo "You can now run the pipeline with:"
echo ""
echo -e "${YELLOW}# Multi-modal test (CRISPR + GEX)${NC}"
echo "nextflow run main.nf \\"
echo "  --fastq_dir data/10x_test_data \\"
echo "  --reference data/references/refdata-gex-GRCh38-2024-A \\"
echo "  --feature_ref data/10x_test_data/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv \\"
echo "  --outdir results/test_run"
echo ""
echo -e "${YELLOW}# BCL conversion test${NC}"
echo "nextflow run main.nf \\"
echo "  --bcl_dir data/miseq-atac \\"
echo "  --samplesheet data/miseq-atac/SampleSheet.csv \\"
echo "  --outdir results/bcl_test \\"
echo "  --skip_alignment --skip_analysis"
echo ""
