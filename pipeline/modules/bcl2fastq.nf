/*
 * BCL2FASTQ MODULE
 * Convert BCL files to FASTQ format (Primary Analysis)
 * Uses Illumina bcl2fastq2 (recommended over deprecated cellranger mkfastq)
 * 
 * Sample Sheet Support:
 * - If provided: Uses sample sheet for demultiplexing samples by index
 * - If not provided: Generates Undetermined_S0_L00X_R1_001.fastq.gz files
 *   (all reads in one file per lane, no demultiplexing)
 */

process BCL2FASTQ {
    tag "${run_id}"
    label 'high_memory'
    publishDir "${params.outdir}/fastq", mode: 'copy'
    
    container 'genexomics/bcl2fastq:1.0'
    
    // Copy BCL directory to work dir to avoid symlink issues in Docker
    stageInMode 'copy'
    
    input:
    tuple val(run_id), path(bcl_dir)
    path samplesheet, stageAs: 'SampleSheet.csv'
    
    output:
    tuple val(run_id), path("${run_id}/outs"), emit: output_dir
    path "${run_id}/outs/fastq_path/**.fastq.gz", emit: fastqs
    path "${run_id}/outs/qc_summary.json", emit: qc_summary
    path "${run_id}/outs/fastq_path", emit: fastq_dir
    path "${run_id}/outs/Reports", emit: reports, optional: true
    path "${run_id}/outs/Stats", emit: stats, optional: true
    
    script:
    // Build sample sheet argument if provided (file is not empty)
    def samplesheet_arg = samplesheet.name != 'NO_SAMPLE_SHEET' ? "--sample-sheet SampleSheet.csv" : ""
    
    """
    # If no sample sheet provided, bcl2fastq will generate Undetermined files
    # using auto-detected settings from RunInfo.xml
    
    bcl2fastq \\
        --runfolder-dir ${bcl_dir} \\
        --output-dir ${run_id}/outs/fastq_path \\
        --reports-dir ${run_id}/outs/Reports \\
        --stats-dir ${run_id}/outs/Stats \\
        ${samplesheet_arg} \\
        --processing-threads ${task.cpus} \\
        --no-lane-splitting \\
        --minimum-trimmed-read-length 35 \\
        --mask-short-adapter-reads 22
    
    # Create QC summary with demultiplexing statistics
    echo '{
      "status": "complete",
      "sample_sheet": "${samplesheet.name}",
      "demultiplexing": "${samplesheet.name != 'NO_SAMPLE_SHEET' ? 'enabled' : 'disabled'}"
    }' > ${run_id}/outs/qc_summary.json
    """
    
    stub:
    """
    mkdir -p ${run_id}/outs/fastq_path
    mkdir -p ${run_id}/outs/Reports
    mkdir -p ${run_id}/outs/Stats
    touch ${run_id}/outs/qc_summary.json
    echo "Stub: bcl2fastq would convert BCL to FASTQ"
    echo "Sample sheet: ${samplesheet.name}"
    """
}
