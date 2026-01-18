/*
 * CELLRANGER_MULTI MODULE
 * Multi-modal analysis for CRISPR + Gene Expression data
 * Processes multiple library types (GEX, CRISPR) together
 */

process CELLRANGER_MULTI {
    tag "${sample_id}"
    label 'high_memory'
    publishDir "${params.outdir}/cellranger_multi", mode: 'copy'
    
    container 'genexomics/cellranger:1.0'
    
    input:
    tuple val(sample_id), path(fastq_dir)
    path reference
    path feature_ref
    path config_csv
    
    output:
    tuple val(sample_id), path("${sample_id}/outs/per_sample_outs/*/sample_filtered_feature_bc_matrix.h5"), emit: matrix
    path "${sample_id}/outs/per_sample_outs/*/web_summary.html", emit: web_summary
    path "${sample_id}/outs/per_sample_outs/*/metrics_summary.csv", emit: summary
    path "${sample_id}/outs/per_sample_outs/*/crispr_analysis", optional: true, emit: crispr_analysis
    path "${sample_id}", emit: full_output
    
    script:
    def feature_ref_path = feature_ref.name != 'NO_FEATURE_REF' ? feature_ref : ''
    def use_provided_config = config_csv.name != 'NO_CONFIG'
    
    // Get absolute paths for Docker container
    def ref_abs_path = reference.toRealPath().toString()
    def feature_ref_abs = feature_ref_path ? feature_ref_path.toRealPath().toString() : ''
    
    """
    # Create multi config CSV with correct paths for this run
    cat > multi_config.csv << EOF
[gene-expression]
reference,${ref_abs_path}
create-bam,true

EOF
    
    # Add feature section if feature ref provided
    if [ -n "${feature_ref_abs}" ]; then
        cat >> multi_config.csv << EOF
[feature]
reference,${feature_ref_abs}

EOF
    fi
    
    # Add libraries section
    cat >> multi_config.csv << 'EOF'
[libraries]
fastq_id,fastqs,feature_types
EOF
    
    # Get absolute path for fastq_dir
    FASTQ_ABS_PATH=\$(cd ${fastq_dir} && pwd)
    
    # Auto-detect library types from subdirectories
    for subdir in ${fastq_dir}/*/; do
        if [ -d "\${subdir}" ]; then
            dirname=\$(basename "\${subdir}")
            dirname_lower=\$(echo "\${dirname}" | tr '[:upper:]' '[:lower:]')
            
            # Get absolute path for this subdirectory
            subdir_abs=\$(cd "\${subdir}" && pwd)
            
            case "\${dirname_lower}" in
                crispr)
                    feature_type="CRISPR Guide Capture"
                    ;;
                gex|gene_expression)
                    feature_type="Gene Expression"
                    ;;
                antibody|adt)
                    feature_type="Antibody Capture"
                    ;;
                multiplexing|multiplexing_capture)
                    feature_type="Multiplexing Capture"
                    ;;
                *)
                    # Default to Gene Expression if cannot determine
                    echo "WARNING: Unknown library type '\${dirname}', defaulting to Gene Expression" >&2
                    feature_type="Gene Expression"
                    ;;
            esac
            
            # Get first FASTQ file to extract sample name
            first_fastq=\$(ls "\${subdir}"/*_R1_*.fastq.gz 2>/dev/null | head -1)
            if [ -n "\${first_fastq}" ]; then
                sample_name=\$(basename "\${first_fastq}" | sed 's/_S[0-9]*_L[0-9]*_R[0-9]*_[0-9]*.fastq.gz//')
                echo "\${sample_name},\${subdir_abs},\${feature_type}" >> multi_config.csv
            fi
        fi
    done
    
    echo "=== Generated Multi Config ===" >&2
    cat multi_config.csv >&2
    echo "=============================" >&2
    
    # Run cellranger multi
    cellranger multi \\
        --id=${sample_id} \\
        --csv=multi_config.csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}
    """
    
    stub:
    """
    mkdir -p ${sample_id}/outs/multi/count
    mkdir -p ${sample_id}/outs/per_sample_outs/${sample_id}
    touch ${sample_id}/outs/multi/count/filtered_feature_bc_matrix.h5
    touch ${sample_id}/outs/multi/count/feature_reference.csv
    touch ${sample_id}/outs/per_sample_outs/${sample_id}/web_summary.html
    echo "Stub: cellranger multi analysis"
    """
}
