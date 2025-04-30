#!/usr/bin/env nextflow

process RSEM_MERGE {

    publishDir "results/rsem_expression", mode: 'copy'
    container 'community.wave.seqera.io/library/r-dplyr_reader:2b13f96b46d7a708'

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