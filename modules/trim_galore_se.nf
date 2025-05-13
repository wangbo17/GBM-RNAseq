#!/usr/bin/env nextflow

process TRIM_GALORE_SE {
    label 'process_medium'
    
    container "containers/trim-galore_0.6.10.sif"
    publishDir "results/trim_galore", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_trimmed.fq.gz", emit: trimmed_reads
    path "*_trimming_report.txt", emit: trimming_reports
    path "*_trimmed_fastqc.{zip,html}", emit: fastqc_reports
    
    script:
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 1
        if (cores < 1) cores = 1
        if (cores > 8) cores = 8
    }

    """
    read_name=\$(basename ${reads})
    target=${sample_id}_raw.fastq.gz

    if [ "\$read_name" != "\$target" ]; then
        ln -s \$read_name \$target
    fi

    trim_galore --fastqc \$target --basename ${sample_id}_val
    """
}
