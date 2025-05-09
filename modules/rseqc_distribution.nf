#!/usr/bin/env nextflow

process RSEQC_DISTRIBUTION {
    label 'process_low'
    
    container 'containers/rseqc_5.0.4_r-base_4.4.3.sif'
    publishDir "results/rseqc", mode: 'copy'

    input:
    path bam_file
    path bed_file

    output:
    path "${bam_file.simpleName}.read_distribution.txt", emit: read_distribution

    script:
    """
    read_distribution.py -i $bam_file -r $bed_file > ${bam_file.simpleName}.read_distribution.txt
    """
}
