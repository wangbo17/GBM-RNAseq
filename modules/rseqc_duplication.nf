#!/usr/bin/env nextflow

process RSEQC_DUPLICATION {
    label 'process_medium'

    container 'containers/rseqc_5.0.4_r-base_4.4.3.sif'
    publishDir "results/rseqc", mode: 'copy'

    input:
    path bam_file

    output:
    path "${bam_file.simpleName}.pos.DupRate.xls", emit: pos_dup_rate
    path "${bam_file.simpleName}.seq.DupRate.xls", emit: seq_dup_rate
    path "${bam_file.simpleName}.DupRate_plot.r", emit: dup_plot_r
    path "${bam_file.simpleName}.DupRate_plot.pdf", emit: dup_plot_pdf

    script:
    """
    read_duplication.py -i $bam_file -o ${bam_file.simpleName}
    """
}
