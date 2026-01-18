#!/usr/bin/env python3
"""
Clustering and Dimensionality Reduction Analysis
================================================
This script performs clustering and visualization of single-cell RNA-seq data.

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
    """Load filtered AnnData object"""
    print(f"Loading data from {input_path}...")
    adata = sc.read_h5ad(input_path)
    print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes")
    return adata


def normalize_and_scale(adata):
    """Normalize and scale data"""
    print("\nNormalizing and scaling data...")
    
    # Normalize to 10,000 counts per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Log-transform
    sc.pp.log1p(adata)
    
    # Scale data (z-score normalization)
    # Use only highly variable genes for scaling
    adata.raw = adata  # Store raw normalized data
    adata = adata[:, adata.var.highly_variable].copy()
    
    sc.pp.scale(adata, max_value=10)
    
    print(f"Scaled {adata.n_vars} highly variable genes")
    
    return adata


def run_pca(adata, n_pcs=50):
    """Perform PCA"""
    print(f"\nRunning PCA with {n_pcs} components...")
    
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
    
    # Plot variance explained (pca_variance_ratio doesn't accept ax parameter)
    sc.pl.pca_variance_ratio(adata, log=True, show=False)
    plt.savefig('pca_variance.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"PCA complete - {n_pcs} components computed")
    
    return adata


def compute_neighborhood_graph(adata, n_pcs=50):
    """Compute neighborhood graph"""
    print(f"\nComputing neighborhood graph using {n_pcs} PCs...")
    
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
    
    print("Neighborhood graph computed")
    
    return adata


def run_umap(adata):
    """Compute UMAP embedding"""
    print("\nComputing UMAP embedding...")
    
    sc.tl.umap(adata)
    
    print("UMAP embedding computed")
    
    return adata


def run_clustering(adata, resolution=0.5):
    """Perform Leiden clustering"""
    print(f"\nPerforming Leiden clustering (resolution={resolution})...")
    
    sc.tl.leiden(adata, resolution=resolution)
    
    n_clusters = len(adata.obs['leiden'].unique())
    print(f"Identified {n_clusters} clusters")
    
    return adata


def plot_clustering_results(adata, sample_id='sample'):
    """Generate clustering visualization plots"""
    print("\nGenerating clustering plots...")
    
    # Combined plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f'Clustering Results - {sample_id}', fontsize=16, y=1.02)
    
    sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, title='Clusters')
    sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[1], show=False, title='Number of Genes')
    
    plt.tight_layout()
    plt.savefig('clustering_results.pdf', bbox_inches='tight')
    plt.close()
    
    # Separate UMAP plot
    fig = sc.pl.umap(adata, color='leiden', show=False, return_fig=True)
    plt.savefig('umap_plot.pdf', bbox_inches='tight')
    plt.close()
    
    print("Clustering plots saved")


def find_marker_genes(adata):
    """Find marker genes for each cluster"""
    print("\nIdentifying marker genes...")
    
    # Find marker genes
    sc.tl.rank_genes_groups(
        adata,
        groupby='leiden',
        method='wilcoxon',
        n_genes=50
    )
    
    # Extract top markers
    markers = []
    for cluster in adata.obs['leiden'].unique():
        cluster_markers = sc.get.rank_genes_groups_df(adata, group=cluster, key='rank_genes_groups')
        cluster_markers['cluster'] = cluster
        markers.append(cluster_markers.head(10))
    
    markers_df = pd.concat(markers)
    markers_df.to_csv('cluster_markers.csv', index=False)
    
    print(f"Saved top 10 marker genes for each cluster")
    
    return adata


def plot_marker_genes(adata):
    """Plot marker gene expression"""
    print("\nPlotting marker genes...")
    
    # Get top 5 markers per cluster
    top_genes = []
    for cluster in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
        markers = sc.get.rank_genes_groups_df(adata, group=cluster, key='rank_genes_groups')
        top_genes.extend(markers.head(2)['names'].tolist())
    
    if len(top_genes) > 0:
        # Limit to 20 genes for readability
        top_genes = top_genes[:20]
        
        fig = sc.pl.dotplot(
            adata,
            var_names=top_genes,
            groupby='leiden',
            show=False,
            return_fig=True
        )
        plt.savefig('marker_genes_dotplot.pdf', bbox_inches='tight')
        plt.close()
        
        print("Marker gene plots saved")


def save_cluster_summary(adata, output_path='cluster_summary.csv'):
    """Save cluster summary statistics"""
    print(f"\nSaving cluster summary to {output_path}...")
    
    cluster_stats = []
    for cluster in sorted(adata.obs['leiden'].unique(), key=lambda x: int(x)):
        cluster_cells = adata.obs['leiden'] == cluster
        stats = {
            'cluster': cluster,
            'n_cells': cluster_cells.sum(),
            'pct_cells': (cluster_cells.sum() / adata.n_obs * 100),
            'median_genes': adata.obs.loc[cluster_cells, 'n_genes_by_counts'].median(),
            'median_counts': adata.obs.loc[cluster_cells, 'total_counts'].median(),
        }
        cluster_stats.append(stats)
    
    pd.DataFrame(cluster_stats).to_csv(output_path, index=False)
    
    print("\nCluster Summary:")
    for stats in cluster_stats:
        print(f"  Cluster {stats['cluster']}: {stats['n_cells']} cells ({stats['pct_cells']:.1f}%)")


def main():
    parser = argparse.ArgumentParser(
        description='Clustering and dimensionality reduction for single-cell RNA-seq'
    )
    parser.add_argument('--input', required=True,
                       help='Input filtered h5ad file')
    parser.add_argument('--output', required=True,
                       help='Output h5ad file with clustering')
    parser.add_argument('--sample-id', default='sample',
                       help='Sample identifier for plots')
    parser.add_argument('--n-pcs', type=int, default=50,
                       help='Number of principal components')
    parser.add_argument('--resolution', type=float, default=0.5,
                       help='Leiden clustering resolution')
    
    args = parser.parse_args()
    
    print("="*70)
    print("Clustering Analysis - GeneXOmics Pipeline")
    print("="*70)
    
    # Load data
    adata = load_data(args.input)
    
    # Normalize and scale
    adata = normalize_and_scale(adata)
    
    # PCA
    adata = run_pca(adata, args.n_pcs)
    
    # Neighborhood graph
    adata = compute_neighborhood_graph(adata, args.n_pcs)
    
    # UMAP
    adata = run_umap(adata)
    
    # Clustering
    adata = run_clustering(adata, args.resolution)
    
    # Visualizations
    plot_clustering_results(adata, args.sample_id)
    
    # Find marker genes
    adata = find_marker_genes(adata)
    
    # Plot markers
    plot_marker_genes(adata)
    
    # Save cluster summary
    save_cluster_summary(adata)
    
    # Save results
    print(f"\nSaving clustered data to {args.output}...")
    adata.write(args.output)
    
    print("\n" + "="*70)
    print("Clustering analysis complete!")
    print("="*70)


if __name__ == '__main__':
    main()
