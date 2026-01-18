/*
 * FASTQC MODULE
 * Quality control analysis of FASTQ files
 */

process FASTQC {
    tag "${sample_id}"
    label 'low_memory'
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    container 'genexomics/fastqc:1.0'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip
    
    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${read1} ${read2}
    """
    
    stub:
    """
    touch ${sample_id}_R1_fastqc.html ${sample_id}_R1_fastqc.zip
    touch ${sample_id}_R2_fastqc.html ${sample_id}_R2_fastqc.zip
    echo "Stub: FastQC quality control analysis"
    """
}
