# Perturb-Seq GeneXOmics Analysis Pipeline

A Nextflow pipeline for processing 10x Genomics CRISPR Perturb-Seq single-cell RNA-seq data, covering primary analysis (BCL to FASTQ conversion) and secondary analysis (alignment, quantification, and clustering).

## Overview

This pipeline processes Perturb-Seq data with automatic data type detection and conditional workflow execution.

### Pipeline Workflow

```
┌─────────────────────────────────────────────────────────────────────┐
│                    INPUT: BCL or FASTQ Files                        │
└────────────────────────────────┬────────────────────────────────────┘
                                 │
                                 ▼
                    ┌────────────────────────┐
                    │  PRIMARY ANALYSIS      │
                    │  • BCL → FASTQ         │ (optional, if BCL input)
                    │  • Quality Control     │
                    └───────────┬────────────┘
                                │
                                ▼
                    ┌────────────────────────┐
                    │  AUTO-DETECT DATA TYPE │
                    └───────────┬────────────┘
                                │
                ┌───────────────┴───────────────┐
                │                               │
                ▼                               ▼
    ┌───────────────────────┐       ┌──────────────────────┐
    │   MULTI-MODAL PATH    │       │  SINGLE-MODAL PATH   │
    │   (CRISPR + GEX)      │       │   (GEX only)         │
    └───────────┬───────────┘       └──────────┬───────────┘
                │                               │
                ▼                               ▼
    ┌───────────────────────┐       ┌──────────────────────┐
    │  Cell Ranger Multi    │       │  Cell Ranger Count   │
    │  • Alignment          │       │  • Alignment         │
    │  • Quantification     │       │  • Quantification    │
    │  • CRISPR assignment  │       └──────────┬───────────┘
    │  • Clustering         │                  │
    │  • UMAP/t-SNE         │                  ▼
    │  • Diff. expression   │       ┌──────────────────────┐
    └───────────┬───────────┘       │  Guide Assignment    │
                │                   └──────────┬───────────┘
                │                              │
                │                              ▼
                │                   ┌──────────────────────┐
                │                   │  QC Filtering        │
                │                   └──────────┬───────────┘
                │                              │
                │                              ▼
                │                   ┌──────────────────────┐
                │                   │  Clustering          │
                │                   │  • PCA, UMAP         │
                │                   │  • Leiden algorithm  │
                │                   └──────────┬───────────┘
                │                              │
                ▼                              ▼
    ┌──────────────────────────────────────────────────────┐
    │           OUTPUTS: Matrices, Reports, Plots          │
    └──────────────────────────────────────────────────────┘
```

**Analysis Stages:**

**Primary Analysis:**

- BCL to FASTQ conversion (bcl2fastq) - optional if starting from BCL
- Quality control (FastQC)

**Secondary Analysis:**

- Read alignment and gene quantification (Cell Ranger count/multi)
- CRISPR guide assignment (multi-modal or count with features)
- Cell quality filtering and QC metrics (single-modal)
- Dimensionality reduction and clustering (single-modal)

## Requirements

- **System**: Linux with x86-64 architecture
- **Resources**: 32 GB RAM (64 GB recommended), 8+ CPU cores, 400 GB disk space
- **Software**: Nextflow ≥21.10.0, Docker ≥20.10.0

## Installation

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
```

### 2. Install Docker

Follow the official guide at <https://docs.docker.com/engine/install/>

```bash
docker --version  # Verify installation
```

### 3. Clone Repository and Build Containers

```bash
git clone <repository-url>
cd pipeline

# Build all Docker containers
docker build -t genexomics/bcl2fastq:1.0 docker/bcl2fastq/
docker build -t genexomics/fastqc:1.0 docker/fastqc/
docker build -t genexomics/cellranger:1.0 docker/cellranger/
docker build -t genexomics/scanpy-qc:1.0 docker/scanpy-qc/
docker build -t genexomics/analysis:1.0 docker/analysis/
```

### 4. Download Test Data and Reference Genome

### Test Dataset Details

**10x Genomics Perturb-Seq Data** (`data/10x_test_data/`)

- **Source**: 10x Genomics 1k CRISPR 5' Gene Expression dataset
- **Type**: Multi-modal (CRISPR + Gene Expression)
- **Size**: ~5 GB
- **Cells**: ~1,000 cells
- **Libraries**: crispr/ and gex/ subdirectories (auto-detected)

**MiSeq BCL Data** (`data/miseq-atac/`)

- **Source**: Illumina MiSeq ATAC-seq run
- **Type**: Raw BCL files with sample sheet
- **Size**: ~500 MB
- **Purpose**: Demonstrates BCL→FASTQ conversion

```bash
# Automated setup (recommended)
bash bin/setup_test_data.sh
```

This script downloads:

- 10x Genomics Perturb-Seq test data (~5 GB)
- MiSeq ATAC BCL test data (~500 MB)
- Cell Ranger reference genome GRCh38-2024-A (~15 GB)

**Total download**: ~20 GB, **estimated time**: 30-60 minutes

## Usage

### Quick Start with Test Data

The pipeline automatically detects whether data is multi-modal (CRISPR+GEX) or single-modal (GEX only).

```bash
# Multi-modal example (CRISPR + Gene Expression)
nextflow run main.nf \
  --fastq_dir data/10x_test_data \
  --reference data/references/refdata-gex-GRCh38-2024-A \
  --feature_ref data/10x_test_data/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv \
  --outdir results/test_run

# Single-modal example (Gene Expression only)
nextflow run main.nf \
  --fastq_dir data/10x_test_data/1k_CRISPR_5p_gemx_fastqs/gex \
  --reference data/references/refdata-gex-GRCh38-2024-A \
  --feature_ref data/10x_test_data/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv \
  --outdir results/test_run
```

**Expected runtime**: 30-60 minutes

**Data Detection:**

- Multi-modal: Subdirectories like `crispr/`, `gex/` → uses Cell Ranger multi
- Single-modal: Flat FASTQ files → uses Cell Ranger count with additional analysis

### Running with own Data

**Starting from BCL files (raw sequencer output):**

```bash
nextflow run main.nf \
  --bcl_dir path/to/bcl/directory \
  --sample_sheet path/to/SampleSheet.csv \
  --reference path/to/reference \
  --feature_ref path/to/feature_reference.csv \
  --outdir results
```

**Starting from FASTQ files:**

```bash
nextflow run main.nf \
  --fastq_dir path/to/fastq/directory \
  --reference path/to/reference \
  --feature_ref path/to/feature_reference.csv \
  --outdir results
```

The pipeline accepts either a run directory (auto-discovers FASTQ subdirectories) or direct path to FASTQ files.

### Key Parameters

| Parameter | Description | Default |
| ----------- | ------------- | --------- |
| `--fastq_dir` / `--bcl_dir` | Input directory | Required |
| `--reference` | Cell Ranger reference genome | Required |
| `--feature_ref` | CRISPR guide reference CSV | Required |
| `--expected_cells` | Expected number of cells | 0 (automatic detection) |
| `--min_genes` | Min genes per cell for QC | 200 |
| `--max_pct_mito` | Max mitochondrial % | 20 |
| `--outdir` | Output directory | results |

## Example Inputs and Outputs

### Input Structure

```
data/
├── 10x_test_data/                         # Run directory (can point here)
│   ├── 1k_CRISPR_5p_gemx_fastqs/          # FASTQ directory (or point here)
│   │   ├── crispr/                        # Multi-modal: CRISPR library
│   │   │   ├── sample_S1_L001_R1_001.fastq.gz
│   │   │   └── sample_S1_L001_R2_001.fastq.gz
│   │   └── gex/                           # Multi-modal: GEX library
│   │       ├── sample_S2_L001_R1_001.fastq.gz
│   │       └── sample_S2_L001_R2_001.fastq.gz
│   └── 1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv
└── references/
    └── refdata-gex-GRCh38-2024-A/
```

**Input Requirements:**

- **FASTQ files**: Paired-end reads following 10x Genomics naming convention (`*_R{1,2}_*.fastq.gz`)
- **Reference genome**: Cell Ranger compatible (downloadable from 10x Genomics)
- **Feature reference**: CSV with format `id,name,read,pattern,sequence,feature_type`

### Output Structure

**For Multi-Modal Data (CRISPR + Gene Expression):**

```
results/
├── fastqc/                             # Quality control reports (~9 MB)
│   ├── {library}_crispr_S1_L001_R1_001_fastqc.html
│   ├── {library}_crispr_S1_L001_R2_001_fastqc.html
│   ├── {library}_gex_S2_L001_R1_001_fastqc.html
│   └── {library}_gex_S2_L001_R2_001_fastqc.html
│
├── cellranger_multi/                   # Complete integrated analysis (3-5 GB)
│   └── {run_id}/
│       ├── outs/
│       │   ├── config.csv                         # Multi config used
│       │   ├── feature_reference.csv              # CRISPR guide reference
│       │   ├── qc_report.html                     # Overall QC report
│       │   ├── qc_library_metrics.csv             # Per-library metrics
│       │   ├── qc_sample_metrics.csv              # Per-sample metrics
│       │   ├── filtered_feature_bc_matrix.h5      # Combined filtered matrix
│       │   ├── raw_feature_bc_matrix.h5           # Combined raw matrix
│       │   ├── raw_cloupe.cloupe                  # Loupe Browser file (all data)
│       │   │
│       │   └── per_sample_outs/
│       │       └── {sample}/
│       │           ├── web_summary.html                        # Interactive QC & analysis
│       │           ├── metrics_summary.csv                     # Key sample metrics
│       │           ├── sample_cloupe.cloupe                    # Loupe Browser (sample)
│       │           ├── sample_filtered_feature_bc_matrix.h5    # Gene + CRISPR matrix
│       │           ├── sample_alignments.bam                   # Aligned reads
│       │           │
│       │           ├── analysis/                               # Clustering & visualization
│       │           │   ├── clustering/
│       │           │   │   └── gene_expression_graphclust/     # Graph-based clusters
│       │           │   ├── pca/                                # PCA projections
│       │           │   ├── umap/                               # UMAP coordinates
│       │           │   ├── tsne/                               # t-SNE coordinates
│       │           │   └── diffexp/                            # Differential expression
│       │           │
│       │           └── crispr_analysis/                        # CRISPR-specific results
│       │               ├── protospacer_calls_per_cell.csv      # Guide assignments
│       │               ├── protospacer_calls_summary.csv       # Assignment statistics
│       │               ├── perturbation_efficiencies_by_feature.csv  # Per-guide effects
│       │               ├── perturbation_efficiencies_by_target.csv   # Per-target effects
│       │               └── cells_per_protospacer.json          # Cell counts per guide
│       │
│       └── [Cell Ranger pipeline metadata files]
│
└── pipeline_info/                      # Execution reports
    ├── report.html                            # Resource usage & metrics (if generated)
    ├── timeline.html                          # Execution timeline
    └── dag.svg                                # Pipeline workflow diagram
```

**Key outputs for multi-modal:**

- `per_sample_outs/{run_id}/sample_filtered_feature_bc_matrix.h5`: Combined gene + CRISPR expression matrix
- `per_sample_outs/{run_id}/web_summary.html`: Interactive QC report with clustering and UMAP
- `per_sample_outs/{run_id}/crispr_analysis/`: CRISPR guide assignments and perturbation effects
- `per_sample_outs/{run_id}/analysis/`: Complete clustering, dimensionality reduction, and differential expression

---

**For Single-Modal Data (Gene Expression Only):**

```
results/
├── fastqc/                             # Quality control reports (4 MB)
│   ├── {sample}_R1_fastqc.html
│   └── {sample}_R2_fastqc.html
│
├── cellranger_count/                   # Cell Ranger outputs (3-5 GB)
│   └── {run_id}/outs/
│       ├── filtered_feature_bc_matrix.h5      # Gene expression matrix
│       ├── web_summary.html                   # QC summary with metrics
│       ├── metrics_summary.csv                # Key quality metrics
│       ├── possorted_genome_bam.bam           # Aligned reads
│       └── cloupe.cloupe                      # Loupe Browser file
│
├── guide_assignment/                   # CRISPR guide assignments (50 MB, if applicable)
│   ├── {sample}_with_guides.h5ad              # Annotated single-cell data
│   ├── {sample}_guide_assignment_summary.csv  # Guide statistics
│   └── {sample}_guide_distribution.pdf        # Visualization plots
│
├── qc_filter/                          # Quality-filtered cells (50 MB)
│   ├── {sample}_filtered.h5ad                 # High-quality cells only
│   ├── {sample}_qc_summary.csv                # QC metrics
│   ├── {sample}_qc_plots.pdf                  # QC visualizations
│   └── {sample}_highly_variable_genes.pdf     # HVG selection plot
│
├── clustering/                         # Final analysis (60 MB)
│   ├── {sample}_clustered.h5ad                # Complete analyzed dataset
│   ├── {sample}_cluster_markers.csv           # Marker genes per cluster
│   ├── {sample}_umap.pdf                      # UMAP visualization
│   └── {sample}_clustering_results.pdf        # Comprehensive results
│
└── pipeline_info/                      # Execution reports
    ├── report.html                            # Resource usage & metrics
    └── timeline.html                          # Execution timeline
```

**Key outputs for single-modal:**

- `filtered_feature_bc_matrix.h5`: Gene expression matrix from Cell Ranger
- `{sample}_clustered.h5ad`: Final analysis with custom clustering parameters
- `{sample}_cluster_markers.csv`: Differentially expressed genes per cluster

## Troubleshooting

**Cell Ranger fails on Apple Silicon (M1/M2 Macs)**  
Cell Ranger requires x86-64 architecture. Use an x86-64 Linux system (e.g., AWS EC2).

**Out of memory errors**  
Cell Ranger requires significant RAM. Use a system with 32+ GB RAM or adjust `--max_memory` parameter.

**Docker not available**  
Ensure Docker daemon is running: `docker ps` should execute without errors.

**Resume interrupted runs**  
Use `-resume` flag to continue from the last successful step:

```bash
nextflow run main.nf <parameters> -resume
```

## Repository Structure

```
pipeline/
├── main.nf                    # Main workflow file
├── nextflow.config            # Configuration and standard values
├── modules/                   # Process definitions
├── docker/                    # Container definitions
├── bin/                       # Helper scripts
│   └── setup_test_data.sh     # Automated test data download
├── conf/
│   └── aws.config             # AWS Batch configuration
└── data/                      # Test data and references
```

## Cloud Deployment

The pipeline is cloud-prepared and can be deployed on AWS using AWS Batch. See [conf/aws.config](conf/aws.config) for:

- AWS Batch executor configuration
- S3 integration for data and work directory
- ECR container registry setup
- Resource specifications per process
- Cost optimization with spot instances
- Retry strategies for reliability

**Quick AWS deployment**:

```bash
# Configure AWS environment
export AWS_REGION=us-east-1
export AWS_BATCH_QUEUE=genexomics-pipeline-queue
export WORK_BUCKET=s3://genexomics-nextflow-work
export ECR_REGISTRY=<account>.dkr.ecr.us-east-1.amazonaws.com

# Run on AWS
nextflow run main.nf -profile aws \
  --fastq_dir s3://genexomics-data/run_001 \
  --reference s3://genexomics-references/refdata-gex-GRCh38-2024-A \
  --outdir s3://genexomics-results/run_001
```

See `conf/aws.config` for detailed setup instructions and cost estimates.

## Docker Containers

All tools are containerized for reproducibility:

- `genexomics/bcl2fastq:1.0`: BCL to FASTQ conversion (bcl2fastq v2.20)
- `genexomics/fastqc:1.0`: Quality control (FastQC v0.12.1)
- `genexomics/cellranger:1.0`: Alignment and quantification (Cell Ranger 10.0.0)
- `genexomics/scanpy-qc:1.0`: Quality filtering (Scanpy 1.9.6)
- `genexomics/analysis:1.0`: Clustering analysis (Scanpy, scikit-learn)

---

**Additional documentation available in repository:**

- `OUTPUT_VALIDATION.md`: Detailed validation report
- `SAMPLE_SHEET_GUIDE.md`: Sample sheet preparation guide
- `conf/aws.config`: AWS configuration profile
