#!/usr/bin/env nextflow

process MULTIQC {
    label 'process_single'
    
    container 'containers/multiqc_1.28.sif'
    publishDir "results/multiqc", mode: 'copy'

    input:
    path '*'
    val output_name

    output:
    path "${output_name}.html", emit: report
    path "${output_name}_data", emit: data

    script:
    """
    multiqc . -n ${output_name}.html
    """
}