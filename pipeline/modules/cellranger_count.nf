/*
 * CELLRANGER_COUNT MODULE
 * Alignment and gene expression quantification (Secondary Analysis)
 * Uses Cell Ranger (industry standard for 10x Genomics data)
 */

process CELLRANGER_COUNT {
    tag "${sample_id}"
    label 'high_memory'
    publishDir "${params.outdir}/cellranger_count", mode: 'copy'
    
    container 'genexomics/cellranger:1.0'
    
    input:
    tuple val(sample_id), path(read1), path(read2), path(fastq_dir)
    path reference
    val expected_cells
    
    output:
    tuple val(sample_id), path("${sample_id}/outs/filtered_feature_bc_matrix.h5"), emit: matrix
    tuple val(sample_id), path("${sample_id}/outs/possorted_genome_bam.bam"), optional: true, emit: bam
    path "${sample_id}/outs/metrics_summary.csv", emit: summary
    path "${sample_id}/outs/web_summary.html", emit: web_summary
    path "${sample_id}/outs/cloupe.cloupe", optional: true, emit: cloupe
    script:
    // Extract sample prefix without lane information for Cell Ranger
    // Cell Ranger expects the prefix (e.g., "sample_name") not "sample_name_S1_L001"
    def sample_prefix = sample_id.replaceAll(/_S\d+.*$/, '')
    
    // Only pass --expect-cells if provided (non-null/non-zero)
    // If skipped, Cell Ranger uses auto-detection (EmptyDrops method) which is generally reliable
    def expect_cells_flag = expected_cells ? "--expect-cells=${expected_cells}" : ""
    
    """
    cellranger count \\
        --id=${sample_id} \\
        --create-bam=true \\
        --transcriptome=${reference} \\
        --fastqs=${fastq_dir} \\
        --sample=${sample_prefix} \\
        ${expect_cells_flag} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}
    """
    
    stub:
    """
    mkdir -p ${sample_id}/outs
    touch ${sample_id}/outs/filtered_feature_bc_matrix.h5
    touch ${sample_id}/outs/metrics_summary.csv
    touch ${sample_id}/outs/web_summary.html
    echo "Stub: cellranger count analysis"
    """
}
