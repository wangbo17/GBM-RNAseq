#!/usr/bin/env nextflow

process RSEQC_BAMSTAT {
    label 'process_low'
    
    container 'containers/rseqc_5.0.4_r-base_4.4.3.sif'
    publishDir "results/rseqc", mode: 'copy'

    input:
    path bam_file

    output:
    path "${bam_file.simpleName}.bam_stat.txt", emit: bamstat_report

    script:
    """
    bam_stat.py -i $bam_file > ${bam_file.simpleName}.bam_stat.txt
    """
}
