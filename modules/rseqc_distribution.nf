#!/usr/bin/env nextflow

process RSEQC_DISTRIBUTION {
    label 'process_low'
    
    container 'oras://community.wave.seqera.io/library/rseqc_r-base:280f3a07db26195c'
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
