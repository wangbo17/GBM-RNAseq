#!/usr/bin/env nextflow

process TRIM_GALORE_PE {
    label 'process_medium'
    
    container "containers/trim-galore_0.6.10.sif"
    publishDir "results/trim_galore", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt", emit: trimming_reports
    path "*_val_1_fastqc.{zip,html}", emit: fastqc_reports_1
    path "*_val_2_fastqc.{zip,html}", emit: fastqc_reports_2

    script:
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 1
        if (cores < 1) cores = 1
        if (cores > 8) cores = 8
    }

    """
    read1_name=\$(basename ${reads[0]})
    read2_name=\$(basename ${reads[1]})

    target1=${sample_id}_raw_1.fastq.gz
    target2=${sample_id}_raw_2.fastq.gz

    if [ "\$read1_name" != "\$target1" ]; then
        ln -s \$read1_name \$target1
    fi

    if [ "\$read2_name" != "\$target2" ]; then
        ln -s \$read2_name \$target2
    fi

    trim_galore --cores $cores --fastqc --paired \$target1 \$target2 --basename ${sample_id}
    """
}
