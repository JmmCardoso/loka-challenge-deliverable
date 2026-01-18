/*
 * GUIDE_ASSIGNMENT MODULE
 * Assigns CRISPR guide RNAs to cells based on feature barcoding data
 * Critical step for Perturb-Seq analysis
 */

process GUIDE_ASSIGNMENT {
    tag "${sample_id}"
    label 'medium_memory'
    publishDir "${params.outdir}/guide_assignment", mode: 'copy'
    
    container 'genexomics/scanpy-qc:1.0'
    
    input:
    tuple val(sample_id), path(matrix_h5)
    path feature_ref
    val assignment_method
    val count_threshold
    
    output:
    tuple val(sample_id), path("${sample_id}_with_guides.h5ad"), emit: h5ad
    path "${sample_id}_guide_assignment_summary.csv", emit: summary
    path "${sample_id}_guide_counts.csv", emit: guide_counts
    path "${sample_id}_guide_distribution.pdf", emit: plots
    
    script:
    """
    #!/usr/bin/env python3
    
    import os
    # Disable Numba caching to avoid permission issues in Docker
    os.environ['NUMBA_CACHE_DIR'] = '/tmp'
    os.environ['NUMBA_DISABLE_JIT'] = '0'
    
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from pathlib import Path
    
    # Configure plotting
    sc.settings.set_figure_params(dpi=100, facecolor='white')
    sns.set_style('whitegrid')
    
    print("Loading data...")
    # Load the Cell Ranger output (includes both gene expression and CRISPR guide counts)
    adata = sc.read_10x_h5('${matrix_h5}')
    
    # Make gene names unique (Cell Ranger may have duplicate gene symbols)
    adata.var_names_make_unique()
    
    # Load feature reference to identify guide features
    feature_df = pd.read_csv('${feature_ref}')
    guide_features = feature_df[feature_df['feature_type'] == 'CRISPR Guide Capture']['id'].tolist()
    
    print(f"Found {len(guide_features)} guide features in reference")
    print(f"Guide features: {guide_features}")
    
    # Separate gene expression and guide counts
    # In 10x Feature Barcoding, guides are in the same matrix as genes
    is_guide = adata.var_names.isin(guide_features)
    
    if is_guide.sum() == 0:
        print("WARNING: No guide features found in the data!")
        print("This might indicate:")
        print("  1. No CRISPR guide reads in the data")
        print("  2. Feature reference doesn't match the data")
        print("  3. Wrong feature_type in reference CSV")
        
        # Create placeholder output for cells with no guides
        adata.obs['guide'] = 'No_Guide_Detected'
        adata.obs['guide_count'] = 0
        adata.obs['n_guides'] = 0
        adata.obs['is_perturbed'] = False
        
    else:
        print(f"Found {is_guide.sum()} guide features in the data")
        
        # Extract guide counts (create a separate view for guides)
        guide_adata = adata[:, is_guide].copy()
        
        # Assignment method: threshold-based or max
        if '${assignment_method}' == 'threshold':
            # Threshold-based: assign guide if count > threshold
            threshold = ${count_threshold}
            
            # For each cell, identify guides above threshold
            guide_assignments = []
            guide_counts_list = []
            n_guides_list = []
            
            for i in range(guide_adata.n_obs):
                cell_guides = guide_adata.X[i, :].toarray().flatten() if hasattr(guide_adata.X, 'toarray') else guide_adata.X[i, :]
                above_threshold = cell_guides > threshold
                
                if above_threshold.sum() == 0:
                    guide_assignments.append('No_Guide')
                    guide_counts_list.append(0)
                    n_guides_list.append(0)
                elif above_threshold.sum() == 1:
                    guide_idx = np.where(above_threshold)[0][0]
                    guide_name = guide_adata.var_names[guide_idx]
                    guide_assignments.append(guide_name)
                    guide_counts_list.append(int(cell_guides[guide_idx]))
                    n_guides_list.append(1)
                else:
                    # Multiple guides detected
                    guide_indices = np.where(above_threshold)[0]
                    guide_names = [guide_adata.var_names[idx] for idx in guide_indices]
                    guide_assignments.append('Multiple_Guides')
                    guide_counts_list.append(int(cell_guides.sum()))
                    n_guides_list.append(len(guide_names))
                    
        else:  # 'max' method
            # Assign the guide with maximum count (if any)
            guide_assignments = []
            guide_counts_list = []
            
            for i in range(guide_adata.n_obs):
                cell_guides = guide_adata.X[i, :].toarray().flatten() if hasattr(guide_adata.X, 'toarray') else guide_adata.X[i, :]
                max_idx = np.argmax(cell_guides)
                max_count = cell_guides[max_idx]
                
                if max_count == 0:
                    guide_assignments.append('No_Guide')
                    guide_counts_list.append(0)
                else:
                    guide_name = guide_adata.var_names[max_idx]
                    guide_assignments.append(guide_name)
                    guide_counts_list.append(int(max_count))
            
            n_guides_list = [1 if g != 'No_Guide' else 0 for g in guide_assignments]
        
        # Add guide assignments to adata.obs
        adata.obs['guide'] = guide_assignments
        adata.obs['guide_count'] = guide_counts_list
        adata.obs['n_guides'] = n_guides_list
        adata.obs['is_perturbed'] = adata.obs['guide'] != 'No_Guide'
    
    # Generate summary statistics
    guide_summary = adata.obs['guide'].value_counts().reset_index()
    guide_summary.columns = ['guide', 'n_cells']
    guide_summary['percentage'] = (guide_summary['n_cells'] / len(adata) * 100).round(2)
    guide_summary.to_csv('${sample_id}_guide_assignment_summary.csv', index=False)
    
    # Save detailed guide counts per cell
    guide_count_df = adata.obs[['guide', 'guide_count', 'n_guides', 'is_perturbed']]
    guide_count_df.to_csv('${sample_id}_guide_counts.csv')
    
    # Generate visualization
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Plot 1: Guide distribution (bar plot)
    ax = axes[0, 0]
    guide_counts = adata.obs['guide'].value_counts().head(20)
    guide_counts.plot(kind='bar', ax=ax, color='steelblue')
    ax.set_title(f'Guide Distribution (Top 20)\\nTotal Cells: {len(adata):,}', fontsize=12, fontweight='bold')
    ax.set_xlabel('Guide', fontsize=10)
    ax.set_ylabel('Number of Cells', fontsize=10)
    ax.tick_params(axis='x', rotation=45, labelsize=8)
    ax.grid(axis='y', alpha=0.3)
    
    # Plot 2: Guide count distribution (histogram)
    ax = axes[0, 1]
    perturbed_cells = adata.obs[adata.obs['is_perturbed']]
    if len(perturbed_cells) > 0:
        ax.hist(perturbed_cells['guide_count'], bins=50, color='coral', edgecolor='black', alpha=0.7)
        ax.set_title(f'Guide UMI Count Distribution\\nPerturbed Cells: {len(perturbed_cells):,}', fontsize=12, fontweight='bold')
        ax.set_xlabel('Guide UMI Count', fontsize=10)
        ax.set_ylabel('Number of Cells', fontsize=10)
        ax.axvline(${count_threshold}, color='red', linestyle='--', linewidth=2, label=f'Threshold ({${count_threshold}})')
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
    else:
        ax.text(0.5, 0.5, 'No perturbed cells detected', ha='center', va='center', fontsize=12)
        ax.set_title('Guide UMI Count Distribution', fontsize=12, fontweight='bold')
    
    # Plot 3: Number of guides per cell
    ax = axes[1, 0]
    n_guides_counts = adata.obs['n_guides'].value_counts().sort_index()
    n_guides_counts.plot(kind='bar', ax=ax, color='mediumseagreen')
    ax.set_title('Number of Guides per Cell', fontsize=12, fontweight='bold')
    ax.set_xlabel('Number of Guides', fontsize=10)
    ax.set_ylabel('Number of Cells', fontsize=10)
    ax.tick_params(axis='x', rotation=0)
    ax.grid(axis='y', alpha=0.3)
    
    # Plot 4: Perturbation rate pie chart
    ax = axes[1, 1]
    perturb_status = adata.obs['is_perturbed'].value_counts()
    
    # Handle case where all cells are unperturbed or all perturbed
    if len(perturb_status) == 1:
        # Only one category present
        if perturb_status.index[0] == True:
            values = [perturb_status.iloc[0], 0]
            labels = [f'Perturbed\\n({perturb_status.iloc[0]:,} cells)', 
                     f'Unperturbed\\n(0 cells)']
        else:
            values = [0, perturb_status.iloc[0]]
            labels = [f'Perturbed\\n(0 cells)', 
                     f'Unperturbed\\n({perturb_status.iloc[0]:,} cells)']
    else:
        # Both categories present
        values = perturb_status.values
        labels = [f'Perturbed\\n({perturb_status.get(True, 0):,} cells)', 
                  f'Unperturbed\\n({perturb_status.get(False, 0):,} cells)']
    
    colors = ['#ff7f0e', '#1f77b4']
    ax.pie(values, labels=labels, autopct='%1.1f%%', colors=colors, 
           startangle=90, textprops={'fontsize': 10})
    ax.set_title('Perturbation Status', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('${sample_id}_guide_distribution.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save annotated AnnData object
    adata.write_h5ad('${sample_id}_with_guides.h5ad')
    
    print(f"\\n{'='*60}")
    print("GUIDE ASSIGNMENT SUMMARY")
    print(f"{'='*60}")
    print(f"Total cells: {len(adata):,}")
    print(f"Perturbed cells: {adata.obs['is_perturbed'].sum():,} ({adata.obs['is_perturbed'].sum()/len(adata)*100:.1f}%)")
    print(f"Unperturbed cells: {(~adata.obs['is_perturbed']).sum():,} ({(~adata.obs['is_perturbed']).sum()/len(adata)*100:.1f}%)")
    print(f"\\nTop 5 guides:")
    print(guide_summary.head(5).to_string(index=False))
    print(f"{'='*60}\\n")
    """
    
    stub:
    """
    touch ${sample_id}_with_guides.h5ad
    touch ${sample_id}_guide_assignment_summary.csv
    touch ${sample_id}_guide_counts.csv
    touch ${sample_id}_guide_distribution.pdf
    echo "Stub: guide assignment analysis"
    """
}
