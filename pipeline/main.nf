#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    GeneXOmics Perturb-Seq Analysis Pipeline
========================================================================================
    Author: João Cardoso
    Date: January 2026
    Description: Automated pipeline for processing 10x Genomics Perturb-Seq data
                 from raw BCL files to clustered single-cell analysis
========================================================================================
*/

// Print pipeline header
def printHeader() {
    log.info """
    ═══════════════════════════════════════════════════════════════════════
     G E N E X O M I C S   P E R T U R B - S E Q   P I P E L I N E
    ═══════════════════════════════════════════════════════════════════════
     Version       : ${workflow.manifest.version}
     Run Name      : ${workflow.runName}
     Profile       : ${workflow.profile}
     
     Input (BCL)   : ${params.bcl_dir ?: params.fastq_dir}
     Output Dir    : ${params.outdir}
     Reference     : ${params.reference}
     Sample Sheet  : ${params.samplesheet ?: 'Not provided'}
     
     Expected Cells: ${params.expected_cells}
     Max Resources : ${params.max_cpus} CPUs, ${params.max_memory} Memory
    ═══════════════════════════════════════════════════════════════════════
    """.stripIndent()
}

// Import modules
include { BCL2FASTQ } from './modules/bcl2fastq'
include { FASTQC } from './modules/fastqc'
include { CELLRANGER_COUNT } from './modules/cellranger_count'
include { CELLRANGER_MULTI } from './modules/cellranger_multi'
include { GUIDE_ASSIGNMENT } from './modules/guide_assignment'
include { QC_FILTER } from './modules/qc_filter'
include { CLUSTERING } from './modules/clustering'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    printHeader()
    
    // Validate parameters
    if (!params.bcl_dir && !params.fastq_dir) {
        error "ERROR: Must provide either --bcl_dir or --fastq_dir"
    }
    
    // Reference is only required if running alignment
    if (!params.reference && !params.skip_alignment) {
        error "ERROR: Must provide --reference (path to reference genome) when running alignment"
    }
    
    // Initialize variables at workflow scope
    def use_multi = false
    def fastq_dir = null
    
    //
    // CHANNEL: Input data
    //
    if (params.bcl_dir) {
        // Start from BCL files (raw sequencer output)
        log.info "Starting from BCL files: ${params.bcl_dir}"
        
        bcl_ch = Channel
            .fromPath(params.bcl_dir, type: 'dir', checkIfExists: true)
            .map { dir -> 
                def run_id = dir.name
                return tuple(run_id, dir)
            }
        
        // Handle optional sample sheet
        // If not provided, bcl2fastq will generate Undetermined files (no demultiplexing)
        if (params.samplesheet) {
            log.info "Using sample sheet: ${params.samplesheet}"
            samplesheet_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
        } else {
            log.warn "⚠️  WARNING: No sample sheet provided"
            log.warn "   bcl2fastq will generate Undetermined_S0 FASTQ files (no demultiplexing)"
            log.warn "   To demultiplex samples, provide a sample sheet with --samplesheet"
            // Create a dummy file channel that signals "no sample sheet"
            samplesheet_ch = Channel.value(file('NO_SAMPLE_SHEET'))
        }
        
        //
        // MODULE: Convert BCL to FASTQ (Primary Analysis)
        //
        BCL2FASTQ(
            bcl_ch,
            samplesheet_ch
        )
        
        // For BCL2FASTQ only run, we just validate the output was created
        // TODO: Parse FASTQ files for downstream processing when needed
        fastq_ch = BCL2FASTQ.out.output_dir
        
    } else {
        // Start from FASTQ files (skip primary analysis)
        log.info "Starting from FASTQ files: ${params.fastq_dir}"
        
        def input_dir = file(params.fastq_dir)
        
        // Validate input directory exists
        if (!input_dir.exists()) {
            error "ERROR: Input directory does not exist: ${params.fastq_dir}\n" +
                  "       Check the path and try again. Use relative paths (e.g., data/10x_test_data) if running locally."
        }
        
        if (!input_dir.isDirectory()) {
            error "ERROR: Input path is not a directory: ${params.fastq_dir}"
        }
        
        // Auto-discover FASTQ directory structure
        // Support both:
        //   1. Direct FASTQ dir: path/to/fastqs/ (contains *.fastq.gz)
        //   2. Run dir: path/to/run_id/ (contains subdirs with FASTQs)
        
        def fastq_files = file("${params.fastq_dir}/*_R{1,2}_*.fastq.gz")
        def has_fastqs_here = fastq_files.size() > 0
        
        fastq_dir = input_dir  // Assign to workflow-level variable
        def run_id = input_dir.name
        
        if (!has_fastqs_here) {
            log.info "📂 No FASTQs in root directory, searching subdirectories..."
            
            // Search one level down for directories containing FASTQs
            def subdirs = input_dir.listFiles()?.findAll { it.isDirectory() } ?: []
            
            if (subdirs.size() == 0) {
                error "ERROR: No subdirectories found in ${params.fastq_dir}\n" +
                      "       Expected either:\n" +
                      "       1. FASTQ files directly in this directory (*_R1_*.fastq.gz, *_R2_*.fastq.gz), OR\n" +
                      "       2. Subdirectories containing FASTQ files"
            }
            
            log.info "   Scanning ${subdirs.size()} subdirectories..."
            subdirs.each { subdir ->
                def count = file("${subdir}/*_R{1,2}_*.fastq.gz").size()
                log.info "     - ${subdir.name}: ${count} FASTQ files"
            }
            
            def fastq_subdir = subdirs.find { subdir ->
                def subdir_fastqs = file("${subdir}/*_R{1,2}_*.fastq.gz")
                if (subdir_fastqs.size() > 0) {
                    return true
                }
                // Or check if it has subdirectories (crispr/, gex/)
                def nested_dirs = subdir.listFiles()?.findAll { it.isDirectory() } ?: []
                return nested_dirs.any { nested -> 
                    file("${nested}/*_R{1,2}_*.fastq.gz").size() > 0 
                }
            }
            
            if (fastq_subdir) {
                log.info "   ✓ Found FASTQ directory: ${fastq_subdir.name}"
                fastq_dir = fastq_subdir
                run_id = input_dir.name  // Keep the parent (run) directory name as run_id
            } else {
                def dir_list = subdirs.collect { it.name }.join(', ')
                error "ERROR: No FASTQ files (*_R{1,2}_*.fastq.gz) found in:\n" +
                      "       - ${params.fastq_dir}\n" +
                      "       - Subdirectories: ${dir_list}\n" +
                      "       Please check:\n" +
                      "       1. FASTQ files follow 10x naming: SampleName_S1_L001_R1_001.fastq.gz\n" +
                      "       2. Files are in the correct location\n" +
                      "       3. You have read permissions"
            }
        }
        
        // Detect if multi-modal data (subdirectories like crispr/, gex/ or multi_config provided)
        def subdirs = fastq_dir.listFiles()?.findAll { it.isDirectory() } ?: []
        def has_multimodal_subdirs = subdirs.any { subdir ->
            def name_lower = subdir.name.toLowerCase()
            def is_known_type = name_lower in ['crispr', 'gex', 'antibody', 'adt', 'gene_expression']
            def has_fastqs = file("${subdir}/*_R{1,2}_*.fastq.gz").size() > 0
            return is_known_type && has_fastqs
        }
        use_multi = params.multi_config != null || has_multimodal_subdirs
        
        if (use_multi) {
            log.info "🔬 Detected multi-modal data (CRISPR + GEX or similar)"
            log.info "   → Using cellranger multi for integrated analysis"
            log.info "   → Run ID: ${run_id}"
            
            // Show detected libraries
            subdirs.findAll { subdir ->
                file("${subdir}/*_R{1,2}_*.fastq.gz").size() > 0
            }.each { subdir ->
                def count = file("${subdir}/*_R{1,2}_*.fastq.gz").size()
                log.info "   → Library: ${subdir.name} (${count} files)"
            }
            
            // For multi, we need the FASTQ directory containing crispr/, gex/, etc.
            fastq_ch = Channel.of(tuple(run_id, fastq_dir))
        } else {
            log.info "📊 Detected single-modal data (GEX only)"
            log.info "   → Using cellranger count for individual samples"
            log.info "   → Run ID: ${run_id}"
            
            def sample_pairs = Channel.fromFilePairs("${fastq_dir}/*_R{1,2}_*.fastq.gz", flat: true)
            def pair_count = file("${fastq_dir}/*_R{1,2}_*.fastq.gz").size() / 2
            log.info "   → Found ${pair_count} sample(s)"
            
            fastq_ch = sample_pairs.map { sample_id, r1, r2 ->
                // Include FASTQ directory for Cell Ranger
                return tuple(sample_id, r1, r2, fastq_dir)
            }
        }
    }
    
    //
    // MODULE: Quality control on FASTQ files
    //
    if (!params.skip_fastqc) {
        if (use_multi) {
            // For multi-modal: create separate channels for each library type
            log.info "📊 Running FastQC on all libraries..."
            
            // Capture fastq_dir value for use in closure
            def fastq_dir_path = fastq_dir.toString()
            def multimodal_subdirs = []
            
            fastq_dir.listFiles()?.each { subdir ->
                if (subdir.isDirectory()) {
                    def fastqs = file("${subdir}/*_R{1,2}_*.fastq.gz")
                    if (fastqs.size() > 0) {
                        multimodal_subdirs.add(subdir.toString())
                    }
                }
            }
            
            fastqc_ch = Channel.fromList(multimodal_subdirs)
                .flatMap { subdir_path ->
                    def pairs = []
                    def r1_files = file("${subdir_path}/*_R1_*.fastq.gz").sort()
                    def r2_files = file("${subdir_path}/*_R2_*.fastq.gz").sort()
                    
                    r1_files.eachWithIndex { r1, idx ->
                        if (idx < r2_files.size()) {
                            def r2 = r2_files[idx]
                            def sample_id = r1.name.replaceAll(/_S\d+_L\d+_R1_\d+\.fastq\.gz$/, '')
                            pairs.add(tuple(sample_id, r1, r2))
                        }
                    }
                    return pairs
                }
            
            FASTQC(fastqc_ch)
        } else {
            // For single-modal: map to remove fastq_dir for FastQC (it only needs sample_id, R1, R2)
            fastqc_ch = fastq_ch.map { sample_id, r1, r2, fastq_dir_param -> tuple(sample_id, r1, r2) }
            FASTQC(fastqc_ch)
        }
    }
    
    //
    // MODULE: Alignment and gene counting (Secondary Analysis)
    //
    if (!params.skip_alignment) {
        reference_ch = Channel.fromPath(params.reference, type: 'dir', checkIfExists: true)
        
        // Route to appropriate Cell Ranger pipeline
        if (use_multi) {
            log.info "📦 Running Cell Ranger Multi pipeline..."
            
            // Prepare inputs for multi
            def multi_config_file = params.multi_config ? 
                Channel.fromPath(params.multi_config, checkIfExists: true) : 
                Channel.value(file('NO_CONFIG'))
            
            def feature_ref_file = params.feature_ref ?
                Channel.fromPath(params.feature_ref, checkIfExists: true) :
                Channel.value(file('NO_FEATURE_REF'))
            
            CELLRANGER_MULTI(
                fastq_ch,
                reference_ch,
                feature_ref_file,
                multi_config_file
            )
            
            log.info "ℹ️  Multi-modal pipeline: Cell Ranger multi provides complete analysis"
            log.info "   ✓ Clustering, UMAP, and differential expression already included"
            log.info "   ✓ CRISPR guide assignment and perturbation analysis completed"
            log.info "   → Skipping redundant QC_FILTER and CLUSTERING modules"
            log.info "   → Final results: ${params.outdir}/cellranger_multi/*/outs/web_summary.html"
            
            // Multi output is already complete - no further analysis needed
            analysis_input = null
            
        } else {
            log.info "📊 Running Cell Ranger Count pipeline..."
            log.info "   → Single-modal data: Running additional QC and clustering analysis"
            
            CELLRANGER_COUNT(
                fastq_ch,
                reference_ch,
                params.expected_cells
            )
            
            //
            // MODULE: Guide Assignment (Perturb-Seq specific)
            //
            if (!params.skip_guide_assignment && params.feature_ref) {
                feature_ref_ch = Channel.fromPath(params.feature_ref, checkIfExists: true)
                
                GUIDE_ASSIGNMENT(
                    CELLRANGER_COUNT.out.matrix,
                    feature_ref_ch,
                    params.guide_assignment_method,
                    params.guide_count_threshold
                )
                
                // Use guide-annotated data for downstream analysis
                analysis_input = GUIDE_ASSIGNMENT.out.h5ad
            } else {
                // Skip guide assignment - use CellRanger output directly
                analysis_input = CELLRANGER_COUNT.out.matrix
            }
        }
        
        //
        // MODULE: Cell quality filtering and QC metrics
        // Only run for single-modal data (Cell Ranger count)
        // Multi-modal data already has complete analysis from Cell Ranger multi
        //
        if (!params.skip_analysis && analysis_input != null) {
            log.info "🔬 Running additional scanpy-based analysis..."
            
            QC_FILTER(
                analysis_input,
                params.min_genes,
                params.min_cells,
                params.max_mito_pct
            )
            
            //
            // MODULE: Clustering and dimensionality reduction
            // Includes: Normalization, PCA, UMAP, Leiden clustering
            //
            CLUSTERING(
                QC_FILTER.out.filtered_h5ad,
                params.n_pcs,
                params.resolution
            )
        } else if (analysis_input == null && !use_multi) {
            log.info "ℹ️  Skipping QC and clustering analysis (--skip_analysis flag set)"
        }
    }
}

/*
========================================================================================
    WORKFLOW INTROSPECTION
========================================================================================
*/

workflow.onComplete {
    def status = workflow.success ? "SUCCESS" : "FAILED"
    def color = workflow.success ? "\033[0;32m" : "\033[0;31m"
    def reset = "\033[0m"
    
    log.info """
    ═══════════════════════════════════════════════════════════════════════
     Pipeline Execution Summary
    ═══════════════════════════════════════════════════════════════════════
     Status        : ${color}${status}${reset}
     Completed at  : ${workflow.complete}
     Duration      : ${workflow.duration}
     Success       : ${workflow.success}
     Exit status   : ${workflow.exitStatus}
     Error report  : ${workflow.errorReport ?: 'No errors'}
     
     Results       : ${params.outdir}
     Work dir      : ${workflow.workDir}
     
     Reports:
       - Timeline  : ${params.outdir}/pipeline_info/timeline.html
       - Report    : ${params.outdir}/pipeline_info/report.html
       - Trace     : ${params.outdir}/pipeline_info/trace.txt
    ═══════════════════════════════════════════════════════════════════════
    """.stripIndent()
}

workflow.onError {
    log.error """
    ═══════════════════════════════════════════════════════════════════════
     Pipeline Execution Error
    ═══════════════════════════════════════════════════════════════════════
     Error message : ${workflow.errorMessage}
     Error report  : ${workflow.errorReport}
    ═══════════════════════════════════════════════════════════════════════
    """.stripIndent()
}
