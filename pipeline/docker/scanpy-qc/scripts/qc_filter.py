#!/usr/bin/env python3
"""
Cell Quality Filtering and QC Analysis
=======================================
This script performs quality control filtering on single-cell RNA-seq data.

Author: João Cardoso
Date: January 2026
"""

import os
# Disable Numba caching to avoid permission issues in Docker
os.environ['NUMBA_CACHE_DIR'] = '/tmp'
os.environ['NUMBA_DISABLE_JIT'] = '0'

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import sys
from pathlib import Path

# Set plotting parameters
sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=False)
sns.set_style("whitegrid")


def load_data(input_path):
    """Load data from various input formats"""
    input_path = Path(input_path)
    
    print(f"Loading data from {input_path}...")
    
    if input_path.suffix == '.h5':
        # Load from Cell Ranger h5 format
        adata = sc.read_10x_h5(input_path)
    elif input_path.suffix == '.h5ad':
        # Load from AnnData format
        adata = sc.read_h5ad(input_path)
    elif input_path.is_dir():
        # Load from MEX format directory
        adata = sc.read_10x_mtx(input_path)
    else:
        raise ValueError(f"Unsupported input format: {input_path}")
    
    print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes")
    return adata


def calculate_qc_metrics(adata):
    """Calculate QC metrics for cells and genes"""
    print("Calculating QC metrics...")
    
    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    # Identify ribosomal genes
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    
    # Identify hemoglobin genes (potential contamination)
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt', 'ribo', 'hb'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    return adata


def plot_qc_metrics(adata, sample_id='sample'):
    """Generate QC plots before filtering"""
    print("Generating QC plots...")
    
    # Create figure with subplots
    fig = plt.figure(figsize=(18, 10))
    fig.suptitle(f'QC Metrics - {sample_id}', fontsize=16, y=0.98)
    
    # Violin plots in top row (3 separate plots)
    ax1 = plt.subplot(2, 3, 1)
    sc.pl.violin(adata, 'n_genes_by_counts', ax=ax1, show=False)
    ax1.set_title('Genes per Cell')
    
    ax2 = plt.subplot(2, 3, 2)
    sc.pl.violin(adata, 'total_counts', ax=ax2, show=False)
    ax2.set_title('Total Counts')
    
    ax3 = plt.subplot(2, 3, 3)
    sc.pl.violin(adata, 'pct_counts_mt', ax=ax3, show=False)
    ax3.set_title('Mitochondrial %')
    
    # Scatter plots in bottom row
    ax4 = plt.subplot(2, 3, 4)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=ax4, show=False)
    
    ax5 = plt.subplot(2, 3, 5)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax5, show=False)
    
    ax6 = plt.subplot(2, 3, 6)
    sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt', ax=ax6, show=False)
    
    plt.tight_layout()
    plt.savefig('qc_violin_plots.pdf', bbox_inches='tight')
    plt.close()
    
    print("QC plots saved")


def filter_cells(adata, min_genes=200, min_cells=3, max_mito_pct=20):
    """Filter low-quality cells and genes"""
    print(f"\nFiltering cells and genes...")
    print(f"  Min genes per cell: {min_genes}")
    print(f"  Min cells per gene: {min_cells}")
    print(f"  Max mitochondrial %: {max_mito_pct}")
    
    # Store initial counts
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Filter cells based on QC metrics
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs.pct_counts_mt < max_mito_pct, :].copy()
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars
    
    print(f"\nFiltering results:")
    print(f"  Cells: {n_cells_before} → {n_cells_after} ({n_cells_after/n_cells_before*100:.1f}% retained)")
    print(f"  Genes: {n_genes_before} → {n_genes_after} ({n_genes_after/n_genes_before*100:.1f}% retained)")
    
    return adata


def identify_highly_variable_genes(adata):
    """Identify highly variable genes"""
    print("\nIdentifying highly variable genes...")
    
    # Normalize and log-transform for HVG detection
    adata_hvg = adata.copy()
    sc.pp.normalize_total(adata_hvg, target_sum=1e4)
    sc.pp.log1p(adata_hvg)
    
    # Identify highly variable genes
    sc.pp.highly_variable_genes(
        adata_hvg,
        n_top_genes=2000,
        flavor='seurat_v3',
        batch_key=None
    )
    
    # Copy HVG information to original adata
    # Note: seurat_v3 creates different columns than seurat flavor
    adata.var['highly_variable'] = adata_hvg.var['highly_variable']
    adata.var['means'] = adata_hvg.var['means']
    
    # seurat_v3 creates 'variances' and 'variances_norm' instead of 'dispersions'
    if 'variances' in adata_hvg.var.columns:
        adata.var['variances'] = adata_hvg.var['variances']
        adata.var['variances_norm'] = adata_hvg.var['variances_norm']
    elif 'dispersions' in adata_hvg.var.columns:
        adata.var['dispersions'] = adata_hvg.var['dispersions']
        adata.var['dispersions_norm'] = adata_hvg.var['dispersions_norm']
    
    n_hvg = adata.var['highly_variable'].sum()
    print(f"Found {n_hvg} highly variable genes")
    
    # Plot highly variable genes
    fig = sc.pl.highly_variable_genes(adata_hvg, show=False)
    plt.savefig('highly_variable_genes.pdf', bbox_inches='tight')
    plt.close()
    
    return adata


def save_qc_summary(adata, output_path='qc_summary.csv'):
    """Save QC summary statistics"""
    print(f"\nSaving QC summary to {output_path}...")
    
    summary = {
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'n_highly_variable_genes': adata.var['highly_variable'].sum(),
        'median_genes_per_cell': adata.obs['n_genes_by_counts'].median(),
        'median_counts_per_cell': adata.obs['total_counts'].median(),
        'median_mito_pct': adata.obs['pct_counts_mt'].median(),
        'mean_genes_per_cell': adata.obs['n_genes_by_counts'].mean(),
        'mean_counts_per_cell': adata.obs['total_counts'].mean(),
        'mean_mito_pct': adata.obs['pct_counts_mt'].mean()
    }
    
    pd.DataFrame([summary]).to_csv(output_path, index=False)
    
    # Print summary
    print("\nQC Summary:")
    for key, value in summary.items():
        print(f"  {key}: {value:.2f}" if isinstance(value, float) else f"  {key}: {value}")


def main():
    parser = argparse.ArgumentParser(
        description='Cell quality filtering and QC for single-cell RNA-seq data'
    )
    parser.add_argument('--input', required=True, 
                       help='Input file (h5, h5ad, or MEX directory)')
    parser.add_argument('--output', required=True,
                       help='Output h5ad file')
    parser.add_argument('--sample-id', default='sample',
                       help='Sample identifier for plots')
    parser.add_argument('--min-genes', type=int, default=200,
                       help='Minimum number of genes per cell')
    parser.add_argument('--min-cells', type=int, default=3,
                       help='Minimum number of cells per gene')
    parser.add_argument('--max-mito', type=float, default=20,
                       help='Maximum mitochondrial percentage')
    
    args = parser.parse_args()
    
    print("="*70)
    print("Cell Quality Filtering - GeneXOmics Pipeline")
    print("="*70)
    
    # Load data
    adata = load_data(args.input)
    
    # Calculate QC metrics
    adata = calculate_qc_metrics(adata)
    
    # Plot QC metrics before filtering
    plot_qc_metrics(adata, args.sample_id)
    
    # Filter cells and genes
    adata = filter_cells(adata, args.min_genes, args.min_cells, args.max_mito)
    
    # Identify highly variable genes
    adata = identify_highly_variable_genes(adata)
    
    # Save QC summary
    save_qc_summary(adata)
    
    # Save filtered data
    print(f"\nSaving filtered data to {args.output}...")
    adata.write(args.output)
    
    print("\n" + "="*70)
    print("QC filtering complete!")
    print("="*70)


if __name__ == '__main__':
    main()
