/*
 * CLUSTERING MODULE
 * Dimensionality reduction and clustering analysis
 */

process CLUSTERING {
    tag "${sample_id}"
    label 'medium_memory'
    publishDir "${params.outdir}/clustering", mode: 'copy'
    
    container 'genexomics/analysis:1.0'
    
    input:
    tuple val(sample_id), path(filtered_h5ad)
    val n_pcs
    val resolution
    
    output:
    tuple val(sample_id), path("${sample_id}_clustered.h5ad"), emit: clustered_h5ad
    path "${sample_id}_clustering_results.pdf", emit: cluster_plots
    path "${sample_id}_pca_variance.pdf", emit: pca_plot
    path "${sample_id}_umap.pdf", emit: umap_plot
    path "${sample_id}_cluster_markers.csv", emit: markers
    
    script:
    """
    python /scripts/clustering.py \\
        --input ${filtered_h5ad} \\
        --output ${sample_id}_clustered.h5ad \\
        --sample-id ${sample_id} \\
        --n-pcs ${n_pcs} \\
        --resolution ${resolution}
    
    # Rename output files with sample ID
    mv clustering_results.pdf ${sample_id}_clustering_results.pdf
    mv pca_variance.pdf ${sample_id}_pca_variance.pdf
    mv umap_plot.pdf ${sample_id}_umap.pdf
    mv cluster_markers.csv ${sample_id}_cluster_markers.csv
    """
    
    stub:
    """
    touch ${sample_id}_clustered.h5ad
    touch ${sample_id}_clustering_results.pdf
    touch ${sample_id}_pca_variance.pdf
    touch ${sample_id}_umap.pdf
    touch ${sample_id}_cluster_markers.csv
    echo "Stub: Clustering analysis"
    """
}
