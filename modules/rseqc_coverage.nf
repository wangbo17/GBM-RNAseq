#!/usr/bin/env nextflow

process RSEQC_COVERAGE {
    label 'process_low'

    container 'containers/rseqc_samtools_r-base.sif'
    publishDir "results/rseqc", mode: 'copy'

    input:
    path bam_file
    path bed_file

    output:
    path "${bam_file.simpleName}.geneBodyCoverage.txt", emit: coverage_report
    path "${bam_file.simpleName}.geneBodyCoverage.r", emit: coverage_rscript
    path "${bam_file.simpleName}.geneBodyCoverage.curves.pdf", emit: coverage_plot

    script:
    """
    samtools index $bam_file
    geneBody_coverage.py -r $bed_file -i $bam_file -o ${bam_file.simpleName}
    """
}
