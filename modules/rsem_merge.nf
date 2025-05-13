#!/usr/bin/env nextflow

process RSEM_MERGE {
    label 'process_low'
    
    container 'containers/r-readr_2.1.5_r-dplyr_1.1.4.sif'
    publishDir "results/rsem_expression", mode: 'copy'

    input:
    path gene_results
    path isoform_results

    output:
    path "gene_expected_count_matrix.tsv"
    path "gene_TPM_matrix.tsv"
    path "gene_FPKM_matrix.tsv"
    path "isoform_expected_count_matrix.tsv"
    path "isoform_TPM_matrix.tsv"
    path "isoform_FPKM_matrix.tsv"
    
    script:
    """
    merge_rsem_matrices.R \\
        --genes ${gene_results.join(' ')} \\
        --isoforms ${isoform_results.join(' ')}
    """
}
