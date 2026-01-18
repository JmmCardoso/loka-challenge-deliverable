/*
 * QC FILTER MODULE
 * Cell quality filtering and QC metrics
 */

process QC_FILTER {
    tag "${sample_id}"
    label 'medium_memory'
    publishDir "${params.outdir}/qc_filter", mode: 'copy'
    
    container 'genexomics/scanpy-qc:1.0'
    
    input:
    tuple val(sample_id), path(matrix_h5)
    val min_genes
    val min_cells
    val max_mito_pct
    
    output:
    tuple val(sample_id), path("${sample_id}_filtered.h5ad"), emit: filtered_h5ad
    path "${sample_id}_qc_summary.csv", emit: qc_summary
    path "${sample_id}_qc_plots.pdf", emit: qc_plots
    path "${sample_id}_highly_variable_genes.pdf", emit: hvg_plot
    
    script:
    """
    python /scripts/qc_filter.py \\
        --input ${matrix_h5} \\
        --output ${sample_id}_filtered.h5ad \\
        --sample-id ${sample_id} \\
        --min-genes ${min_genes} \\
        --min-cells ${min_cells} \\
        --max-mito ${max_mito_pct}
    
    # Rename output plots with sample ID
    mv qc_violin_plots.pdf ${sample_id}_qc_plots.pdf
    mv qc_summary.csv ${sample_id}_qc_summary.csv
    mv highly_variable_genes.pdf ${sample_id}_highly_variable_genes.pdf
    """
    
    stub:
    """
    touch ${sample_id}_filtered.h5ad
    touch ${sample_id}_qc_summary.csv
    touch ${sample_id}_qc_plots.pdf
    touch ${sample_id}_highly_variable_genes.pdf
    echo "Stub: QC filtering"
    """
}
