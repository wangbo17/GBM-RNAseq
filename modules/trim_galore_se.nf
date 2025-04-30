#!/usr/bin/env nextflow

process TRIM_GALORE_SE {

    container "community.wave.seqera.io/library/trim-galore:0.6.10--bc38c9238980c80e"
    publishDir "results/trim_galore", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_trimmed.fq.gz", emit: trimmed_reads
    path "*_trimming_report.txt", emit: trimming_reports
    path "*_trimmed_fastqc.{zip,html}", emit: fastqc_reports

    script:
    """
    read_name=\$(basename ${reads})
    target=${sample_id}_raw.fastq.gz

    if [ "\$read_name" != "\$target" ]; then
        ln -s \$read_name \$target
    fi

    trim_galore --fastqc \$target --basename ${sample_id}_val
    """
}
