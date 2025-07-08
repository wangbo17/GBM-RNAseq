#!/usr/bin/env nextflow

process RSEM_MERGE {
    label 'process_low'
    
    container 'oras://community.wave.seqera.io/library/r-dplyr_r-readr:71c011f9534aeaa9'
    publishDir "results/rsem_expression", mode: 'copy'

    input:
    path gene_results
    path isoform_results

    output:
    path "gene_expected_count_matrix.csv"
    path "gene_TPM_matrix.csv"
    path "gene_FPKM_matrix.csv"
    path "isoform_expected_count_matrix.csv"
    path "isoform_TPM_matrix.csv"
    path "isoform_FPKM_matrix.csv"
    
    script:
    """
    merge_rsem_matrices.R \\
        --genes ${gene_results.join(' ')} \\
        --isoforms ${isoform_results.join(' ')}
    """
}
