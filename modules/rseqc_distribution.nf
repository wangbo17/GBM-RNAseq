#!/usr/bin/env nextflow

process RSEQC_DISTRIBUTION {

    container "community.wave.seqera.io/library/rseqc:5.0.4--8c855f0b915334d2"
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
