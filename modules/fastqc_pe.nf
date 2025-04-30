#!/usr/bin/env nextflow

process FASTQC_PE {

    container "community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29"
    publishDir "results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
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

    fastqc \$target1 \$target2
    """
}