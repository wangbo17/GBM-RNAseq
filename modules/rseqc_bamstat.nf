#!/usr/bin/env nextflow

process RSEQC_BAMSTAT {

    container "community.wave.seqera.io/library/rseqc:5.0.4--8c855f0b915334d2"
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
