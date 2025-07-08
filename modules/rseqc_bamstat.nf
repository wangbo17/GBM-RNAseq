#!/usr/bin/env nextflow

process RSEQC_BAMSTAT {
    label 'process_low'
    
    container 'oras://community.wave.seqera.io/library/rseqc_r-base:280f3a07db26195c'
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
