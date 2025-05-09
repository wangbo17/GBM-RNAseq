#!/usr/bin/env nextflow

process FASTQC_SE {
    label 'process_single'
    
    container 'containers/fastqc_0.12.1.sif'
    publishDir "results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    read_name=\$(basename ${reads})
    target=${sample_id}_raw.fastq.gz

    if [ "\$read_name" != "\$target" ]; then
        ln -s \$read_name \$target
    fi

    fastqc \$target
    """
}